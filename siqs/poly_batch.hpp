#pragma once

// ---------------------------------------------------------------------------
// siqs/poly_batch.hpp — сито одного A-полинома (горячий цикл).
//
// update_roots сдвигает корни при переходе к следующему полиному,
// process_one_A_batch собирает порцию полиномов одного A: просев по байтовым
// логарифмам, отбор выживших, точное разложение кофактора и выдача
// соотношений.
// ---------------------------------------------------------------------------

#include "common.hpp"
#include "relation.hpp"
#include "worker_pool.hpp"
#include "config.hpp"

namespace bignum::siqs_xor {

// Сдвиг корней при переходе к следующему полиному: x += delta (add = true) или
// x -= delta (add = false) по модулю p. Условие `if (inv_A_mod_p[idx])` убрано:
// у выбранных в A простых delta = 0, поэтому сдвиг нулевой в обоих режимах.
// Скалярная ветка на каждом простом стоила дороже самой арифметики.
inline void update_roots(u32 num_primes, const u32 *__restrict prime,
                         const u32 *__restrict delta, u32 *__restrict sol1,
                         u32 *__restrict sol2, bool add)
{
    u32 idx = 1;
#if defined(__SSE2__)
    if (add)
    {
        // sol < p и delta < p, значит сумма < 2p < 2^31: знаковое сравнение
        // корректно, вычитание p делается по маске.
        for (; idx + 4 <= num_primes; idx += 4)
        {
            const __m128i d = _mm_loadu_si128(reinterpret_cast<const __m128i *>(delta + idx));
            const __m128i pp = _mm_loadu_si128(reinterpret_cast<const __m128i *>(prime + idx));
            __m128i a = _mm_add_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol1 + idx)), d);
            __m128i b = _mm_add_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol2 + idx)), d);
            a = _mm_sub_epi32(a, _mm_andnot_si128(_mm_srai_epi32(_mm_sub_epi32(a, pp), 31), pp));
            b = _mm_sub_epi32(b, _mm_andnot_si128(_mm_srai_epi32(_mm_sub_epi32(b, pp), 31), pp));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol1 + idx), a);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol2 + idx), b);
        }
    }
    else
    {
        for (; idx + 4 <= num_primes; idx += 4)
        {
            const __m128i d = _mm_loadu_si128(reinterpret_cast<const __m128i *>(delta + idx));
            const __m128i pp = _mm_loadu_si128(reinterpret_cast<const __m128i *>(prime + idx));
            __m128i a = _mm_sub_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol1 + idx)), d);
            __m128i b = _mm_sub_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol2 + idx)), d);
            a = _mm_add_epi32(a, _mm_and_si128(pp, _mm_srai_epi32(a, 31)));
            b = _mm_add_epi32(b, _mm_and_si128(pp, _mm_srai_epi32(b, 31)));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol1 + idx), a);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol2 + idx), b);
        }
    }
#endif
    if (add)
    {
        for (; idx < num_primes; ++idx)
        {
            const u32 p = prime[idx];
            const u32 d = delta[idx];
            const u32 s1_n = sol1[idx] + d;
            sol1[idx] = (s1_n >= p) ? (s1_n - p) : s1_n;
            const u32 s2_n = sol2[idx] + d;
            sol2[idx] = (s2_n >= p) ? (s2_n - p) : s2_n;
        }
    }
    else
    {
        for (; idx < num_primes; ++idx)
        {
            const u32 p = prime[idx];
            const u32 d = delta[idx];
            sol1[idx] = (sol1[idx] >= d) ? (sol1[idx] - d) : (sol1[idx] + p - d);
            sol2[idx] = (sol2[idx] >= d) ? (sol2[idx] - d) : (sol2[idx] + p - d);
        }
    }
}

// Один батч A-полиномов (полностью однопоточный).
// raw_found/raw_stop_at — общий счётчик сырых соотношений и порог ранней
// остановки: рабочий перестаёт генерировать полиномы, как только батч набрал
// достаточно. nullptr/0 — ограничение выключено.
inline BatchResult process_one_A_batch(const SIQSConfig &cfg, u64 seed,
                                       std::atomic<u64> *raw_found = nullptr,
                                       u64 raw_stop_at = 0, int poly_limit = 0)
{
    BatchResult result;
    SieveStats &st = result.stats;
    auto &found_relations = result.relations;
    auto clock_now = []() { return std::chrono::steady_clock::now(); };
    auto since = [](std::chrono::steady_clock::time_point a,
                    std::chrono::steady_clock::time_point b)
    { return std::chrono::duration<double>(b - a).count(); };
    auto emit_relation = [&](std::array<u64, 3> lp, uint8_t n_lp, const Big &x,
                             int sign, const std::vector<FactorExp> &e)
    {
        found_relations.push_back({lp, n_lp, x, sign, e});
        if (raw_found)
            raw_found->fetch_add(1, std::memory_order_relaxed);
    };
    auto t_phase = clock_now();

    const u32 num_primes = cfg.num_primes;
    const u32 B_size = cfg.B_size;
    const u32 half_B = cfg.half_B;
    const int s = cfg.s;

    std::mt19937_64 rng(seed);
    std::vector<u32> chosen_indices;
    std::vector<u32> pool_copy = cfg.pool_indices;
    for (int i = 0; i < s; ++i)
    {
        std::uniform_int_distribution<size_t> dist(i, pool_copy.size() - 1);
        const size_t pick = dist(rng);
        std::swap(pool_copy[i], pool_copy[pick]);
        chosen_indices.push_back(pool_copy[i]);
    }
    std::sort(chosen_indices.begin(), chosen_indices.end());

    std::vector<uint8_t> is_chosen(num_primes, 0);
    Big A(1);
    for (u32 idx : chosen_indices)
    {
        is_chosen[idx] = 1;
        A = A * Big(cfg.prime_base[idx]);
    }

    std::vector<Big> B_terms(s);
    for (int j = 0; j < s; ++j)
    {
        const u32 idx = chosen_indices[j];
        const u64 q_j = cfg.prime_base[idx];
        const u64 r_j = cfg.root1[idx];
        Big A_div_q = A;
        A_div_q.divrem_u64(q_j);
        const u64 rem_a = A_div_q.mod_u64(q_j);
        const u64 inv = inv_mod64(rem_a, q_j);
        u64 gamma = mul_mod64(r_j, inv, q_j);
        if (gamma > q_j / 2)
            gamma = q_j - gamma;
        B_terms[j] = A_div_q * Big(gamma);
    }

    // B = sum(+-B_terms[j]) — держим модуль и знак отдельно.
    Big B_mag(0);
    bool B_neg = false;
    for (int j = 0; j < s; ++j)
        B_mag.add(B_terms[j]);
    auto b_signed_add = [&](const Big &term, bool term_neg)
    {
        if (B_neg == term_neg)
        {
            B_mag = B_mag + term;
        }
        else if (B_mag >= term)
        {
            B_mag = B_mag - term;
        }
        else
        {
            B_mag = term - B_mag;
            B_neg = term_neg;
        }
        if (B_mag.is_zero())
            B_neg = false;
    };

    std::vector<u32> inv_A_mod_p(num_primes, 0);
    std::vector<std::vector<u32>> B_terms_mod_p(s, std::vector<u32>(num_primes, 0));
    std::vector<u32> sol1(num_primes, 0);
    std::vector<u32> sol2(num_primes, 0);

    // floor(2^32 / p) для быстрого остатка в пробном делении: на 300 битах
    // s_idx % p считается для ~6 тыс. простых на каждого кандидата, а деление
    // на переменный делитель — это 20+ тактов. Умножение с восстановлением
    // остатка стоит ~5 тактов и даёт тот же результат.
    std::vector<u32> p_recip(num_primes, 0);

    for (u32 idx = 1; idx < num_primes; ++idx)
    {
        const u32 p = cfg.prime_base[idx];
        p_recip[idx] = static_cast<u32>((1ULL << 32) / p);
        if (is_chosen[idx])
            continue;
        const u64 a_mod = A.mod_u64(p);
        const u32 inv_A = static_cast<u32>(inv_mod64(a_mod, p));
        inv_A_mod_p[idx] = inv_A;
        for (int j = 0; j < s; ++j)
        {
            const u64 b_mod = B_terms[j].mod_u64(p);
            B_terms_mod_p[j][idx] = static_cast<u32>((2ULL * b_mod * inv_A) % p);
        }

        u64 B_mod = B_mag.mod_u64(p);
        if (B_neg && B_mod != 0)
            B_mod = static_cast<u64>(p) - B_mod;
        const u64 r1 = cfg.root1[idx];
        const u64 r2 = cfg.root2[idx];
        sol1[idx] = static_cast<u32>((((p + r1 - B_mod) % p) * static_cast<u64>(inv_A)
                                      + static_cast<u64>(half_B)) % p);
        sol2[idx] = static_cast<u32>((((p + r2 - B_mod) % p) * static_cast<u64>(inv_A)
                                      + static_cast<u64>(half_B)) % p);
    }

    std::vector<int> signs(s, 1);
    const int polys_per_A = 1 << (s - 1);
    // Ограничение числа полиномов на один A (poly_limit) позволяет набирать
    // соотношения малыми порциями: на небольших числах один полный батч даёт
    // десятки тысяч соотношений, и разреженный Гаусс потом считается впустую.
    const int polys_to_run = (poly_limit <= 0 || poly_limit > polys_per_A)
                                 ? polys_per_A
                                 : poly_limit;

    std::vector<uint8_t> sieve_logs(B_size);
    // Попадания крупных простых в блок: старшие 32 бита — позиция, младшие —
    // индекс простого. Раньше здесь лежали индексы, и второй проход (отбор
    // выживших) заново читал sol1/sol2 по случайным адресам: на 300 битах это
    // десятки тысяч промахов кэша на полином. Теперь проход идёт по плотному
    // массиву пар, а индексация нужна только для редких выживших.
    std::vector<u64> in_block_large1;
    std::vector<u64> in_block_large2;
    in_block_large1.reserve(4096);
    in_block_large2.reserve(4096);

    // Битовая карта выживших позиций и список попаданий крупных простых.
    // Хеш-таблица на полином здесь была лишней: выживших единицы, а проверок —
    // десятки тысяч, и каждая стоила похода в unordered_map. Карта на блок
    // в 128 КБ занимает 16 КБ и целиком лежит в L1.
    std::vector<u64> survivor_bits((B_size + 63) / 64);
    std::vector<u32> survivor_idx;
    // Отобранные пары (позиция, индекс простого) вместо двух параллельных
    // векторов: одна загрузка на проверку в цикле разложения.
    std::vector<u64> hit_pairs;
    std::vector<u32> large_divs;
    survivor_idx.reserve(1024);
    hit_pairs.reserve(64);
    large_divs.reserve(16);
    // Переиспользуемый буфер множителей: раньше вектор создавался на каждого
    // выжившего кандидата (сотни тысяч аллокаций за проход).
    std::vector<FactorExp> exps;
    exps.reserve(64);

    const int base_target_bits = static_cast<int>(cfg.kM.bit_width() / 2
                                                  + std::countr_zero(half_B));
    const int lp2_bits = static_cast<int>(std::bit_width(cfg.lp2_bound));
    const int cof_bits = lp2_bits + static_cast<int>(cfg.lp3_slack_bits);
    const uint8_t target_log = static_cast<uint8_t>(
        std::min(250, std::max(15, base_target_bits - cof_bits + 2)));

    st.t_prep += since(t_phase, clock_now());

    for (int poly_i = 0; poly_i < polys_to_run; ++poly_i)
    {
        // Ранняя остановка. Порция полиномов не должна разрастаться сверх
        // нужного: на малых числах один полином даёт сотни соотношений, и
        // слепое удвоение poly_limit собирало матрицу в разы больше нужной
        // (отсюда 100 бит медленнее 128 в прежних замерах).
        if (poly_i > 0 && raw_found
            && raw_found->load(std::memory_order_relaxed) >= raw_stop_at)
            break;

        if ((poly_i & 15) == 0 && stop_requested())
            return result;

        ++st.polynomials;
        t_phase = clock_now();

        if (poly_i > 0)
        {
            const int gray_curr = poly_i ^ (poly_i >> 1);
            const int gray_prev = (poly_i - 1) ^ ((poly_i - 1) >> 1);
            const int diff = gray_curr ^ gray_prev;
            const int v = static_cast<int>(std::countr_zero(static_cast<uint32_t>(diff))) + 1;

            const int new_sign = -signs[v];
            signs[v] = new_sign;
            const auto &row_v = B_terms_mod_p[v];
            Big two_Bv = B_terms[v];
            two_Bv.shl1();

            if (new_sign == -1)
            {
                b_signed_add(two_Bv, true);
                update_roots(num_primes, cfg.prime_base.data(), row_v.data(), sol1.data(),
                             sol2.data(), true);
            }
            else
            {
                b_signed_add(two_Bv, false);
                update_roots(num_primes, cfg.prime_base.data(), row_v.data(), sol1.data(),
                             sol2.data(), false);
            }
        }

        std::memset(sieve_logs.data(), 0, B_size);
        uint8_t *__restrict s_ptr = sieve_logs.data();
        st.t_prep += since(t_phase, clock_now());
        t_phase = clock_now();

        // 1. Просеивание малых и средних простых.
        for (u32 idx = cfg.small_prime_skip; idx < cfg.split_idx; ++idx)
        {
            if (inv_A_mod_p[idx] == 0)
                continue;
            const u32 p = cfg.prime_base[idx];
            const uint8_t lp = cfg.log_primes[idx];
            const u32 s1 = sol1[idx];
            const u32 s2 = sol2[idx];
            for (u32 s_idx = s1; s_idx < B_size; s_idx += p)
                s_ptr[s_idx] += lp;
            if (s2 != s1)
                for (u32 s_idx = s2; s_idx < B_size; s_idx += p)
                    s_ptr[s_idx] += lp;
        }

        // 2. Просеивание крупных простых (p >= B_size).
        in_block_large1.clear();
        in_block_large2.clear();
        for (u32 idx = cfg.split_idx; idx < num_primes; ++idx)
        {
            const u32 s1 = sol1[idx];
            const u32 s2 = sol2[idx];
            const uint8_t lp = cfg.log_primes[idx];
            if (s1 < B_size)
            {
                s_ptr[s1] += lp;
                in_block_large1.push_back((static_cast<u64>(s1) << 32) | idx);
            }
            if (s2 < B_size)
            {
                s_ptr[s2] += lp;
                in_block_large2.push_back((static_cast<u64>(s2) << 32) | idx);
            }
        }

        st.large_prime_hits += in_block_large1.size() + in_block_large2.size();
        st.t_sieve += since(t_phase, clock_now());
        t_phase = clock_now();

        // 3. Поиск выживших кандидатов.        // Скан блока — самая дорогая часть после просева: на полином приходится
        // B_size побайтовых сравнений. Сравниваем по 16 байт (SSE2), а результат
        // сразу превращаем в индексы через movemask. Байты логарифмов знаковые
        // только формально, поэтому перед cmpgt сдвигаем диапазон на 0x80.
        survivor_idx.clear();
        std::fill(survivor_bits.begin(), survivor_bits.end(), 0ULL);
#if defined(__SSE2__)
        {
            const __m128i threshold =
                _mm_set1_epi8(static_cast<char>((target_log - 1) ^ 0x80));
            const __m128i flip = _mm_set1_epi8(static_cast<char>(0x80));
            u32 s_idx = 0;
            for (; s_idx + 16 <= B_size; s_idx += 16)
            {
                __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s_ptr + s_idx));
                v = _mm_xor_si128(v, flip);
                unsigned mask =
                    static_cast<unsigned>(_mm_movemask_epi8(_mm_cmpgt_epi8(v, threshold)));
                while (mask)
                {
                    survivor_idx.push_back(s_idx + static_cast<u32>(std::countr_zero(mask)));
                    mask &= mask - 1;
                }
            }
            for (; s_idx < B_size; ++s_idx)
                if (s_ptr[s_idx] >= target_log)
                    survivor_idx.push_back(s_idx);
        }
#else
        for (u32 s_idx = 0; s_idx < B_size; ++s_idx)
            if (s_ptr[s_idx] >= target_log)
                survivor_idx.push_back(s_idx);
#endif

        st.survivors += survivor_idx.size();
        if (survivor_idx.empty())
        {
            st.t_scan += since(t_phase, clock_now());
            continue;
        }

        for (u32 s_idx : survivor_idx)
            survivor_bits[s_idx >> 6] |= 1ULL << (s_idx & 63);

        // Крупные простые (p >= B_size) делят Q(x) только там, где их корень
        // попал в блок. Отбираем попадания в выжившие позиции: проверка — один
        // бит из карты в L1 вместо поиска в хеш-таблице на каждое простое.
        // Позиция уже лежит рядом с индексом, поэтому читать sol1/sol2 заново
        // не нужно — проход идёт по плотному массиву.
        hit_pairs.clear();
        for (u64 h : in_block_large1)
        {
            const u32 pos = static_cast<u32>(h >> 32);
            if ((survivor_bits[pos >> 6] >> (pos & 63)) & 1ULL)
                hit_pairs.push_back(h);
        }
        for (u64 h : in_block_large2)
        {
            const u32 pos = static_cast<u32>(h >> 32);
            if ((survivor_bits[pos >> 6] >> (pos & 63)) & 1ULL)
                hit_pairs.push_back(h);
        }

        st.hits_in_survivors += hit_pairs.size();
        st.t_scan += since(t_phase, clock_now());
        t_phase = clock_now();

        // 4. Точное разложение выживших кандидатов.
        for (u32 s_idx : survivor_idx)
        {
            large_divs.clear();
            for (u64 h : hit_pairs)
                if (static_cast<u32>(h >> 32) == s_idx)
                    large_divs.push_back(static_cast<u32>(h));

            const u32 mask = half_B ^ s_idx;
            const i64 y = static_cast<i64>(half_B ^ mask) - static_cast<i64>(half_B);

            const bool y_neg = (y < 0);
            const Big ay = A * Big(static_cast<u64>(y_neg ? -y : y));
            Big X_raw;
            bool x_neg;
            if (y_neg == B_neg)
            {
                X_raw = ay + B_mag;
                x_neg = y_neg;
            }
            else if (ay >= B_mag)
            {
                X_raw = ay - B_mag;
                x_neg = y_neg;
            }
            else
            {
                X_raw = B_mag - ay;
                x_neg = B_neg;
            }

            const Big X2 = X_raw * X_raw;
            if (X2 == cfg.kM)
                continue;

            int sign = 1;
            Big diff;
            if (X2 >= cfg.kM)
            {
                diff = X2 - cfg.kM;
                sign = 1;
            }
            else
            {
                diff = cfg.kM - X2;
                sign = -1;
            }

            Big temp = diff / A;
            exps.clear();
            for (u32 a_idx : chosen_indices)
                exps.push_back({a_idx, 1});

            if ((temp.low_word() & 1) == 0)
            {
                u32 cnt2 = 0;
                while (!temp.is_zero() && (temp.low_word() & 1) == 0)
                {
                    const u32 step = (temp.low_word() != 0)
                                         ? std::countr_zero(temp.low_word())
                                         : 64U;
                    temp.shr(step);
                    cnt2 += step;
                }
                if (cnt2 != 0)
                    exps.push_back({0, cnt2});
            }

            for (u32 idx : large_divs)
            {
                const u32 p = cfg.prime_base[idx];
                if (temp.mod_u64(p) == 0)
                {
                    u32 cnt = 0;
                    while (temp.mod_u64(p) == 0)
                    {
                        temp.divrem_u64(p);
                        ++cnt;
                    }
                    merge_exp(exps, idx, cnt);
                }
            }

            for (u32 idx = 1; idx < cfg.split_idx; ++idx)
            {
                const u32 p = cfg.prime_base[idx];
                // Остаток через предвычисленное floor(2^32 / p) вместо деления.
                const u32 q_est = static_cast<u32>((static_cast<u64>(s_idx) * p_recip[idx]) >> 32);
                u32 rem_s = s_idx - q_est * p;
                if (rem_s >= p)
                    rem_s -= p;
                if (rem_s == sol1[idx] || rem_s == sol2[idx] || is_chosen[idx])
                {
                    if (temp.mod_u64(p) == 0)
                    {
                        u32 cnt = 0;
                        while (temp.mod_u64(p) == 0)
                        {
                            temp.divrem_u64(p);
                            ++cnt;
                        }
                        merge_exp(exps, idx, cnt);
                        if (temp.bit_width() <= 1)
                            break;
                    }
                }
            }

            if (temp.bit_width() > 63)
            {
                ++st.dropped_big_cofactor;
                continue;
            }
            u64 cof = temp.low_word();
            // Без 3LP верхняя граница кофактора прежняя — lp2_bound: иначе
            // лишние кандидаты уходили бы в разложение впустую.
            const u64 cof_gate = cfg.lp3_slack_bits ? cfg.lp3_bound : cfg.lp2_bound;
            if (cof > cof_gate)
            {
                ++st.dropped_big_cofactor;
                continue;
            }

            if (cof > 1 && cof <= cfg.last_prime)
            {
                auto pit = cfg.prime_to_idx.find(cof);
                if (pit != cfg.prime_to_idx.end())
                {
                    merge_exp(exps, pit->second, 1);
                    cof = 1;
                }
                else
                {
                    ++st.dropped_nonsplit;
                    continue;
                }
            }

            Big X_val = X_raw % cfg.M;
            if (x_neg && !X_val.is_zero())
                X_val = cfg.M - X_val;

            if (cof == 1)
            {
                ++st.rel_full;
                emit_relation({{1, 1, 1}}, 0, X_val, sign, exps);
            }
            else if (cof <= cfg.lp_bound)
            {
                ++st.rel_1lp;
                emit_relation({{cof, 1, 1}}, 1, X_val, sign, exps);
            }
            else if (cof >= cfg.min_2lp_comp)
            {
                if (is_prime64(cof))
                {
                    // Одно большое простое за пределами lp_bound: пара не
                    // собирается, соотношение бесполезно.
                    ++st.dropped_big_cofactor;
                }
                else
                {
                    u64 f1 = 0, f2 = 0;
                    if (brent_factor_two_primes(cof, f1, f2) && is_prime64(f1)
                        && is_prime64(f2) && f1 > cfg.last_prime && f1 <= cfg.lp_bound
                        && f2 > cfg.last_prime && f2 <= cfg.lp_bound)
                    {
                        ++st.rel_2lp;
                        emit_relation({{f1, f2, 1}}, 2, X_val, sign, exps);
                    }
                    else
                    {
                        // Кофактор не делится на два допустимых простых —
                        // пробуем три. Пока только считаем: в матрицу такие
                        // соотношения не идут и на раннюю остановку не влияют,
                        // иначе сдвинулась бы траектория A-наборов.
                        u64 g1 = 0, g2 = 0, g3 = 0;
                        if (cfg.lp3_slack_bits
                            && brent_factor_three_primes(cof, cfg.last_prime,
                                                         cfg.lp_bound, g1, g2, g3))
                            ++st.rel_3lp;
                        else
                            ++st.dropped_nonsplit;
                    }
                }
            }
            else
            {
                ++st.dropped_nonsplit;
            }
        }

        st.t_factor += since(t_phase, clock_now());
    }
    return result;
}

} // namespace bignum::siqs_xor
