#pragma once

// ---------------------------------------------------------------------------
// siqs/driver.hpp — автоподбор параметров и главный цикл решета.
//
// auto_tune_siqs / find_best_multiplier / sieve_primes готовят параметры под
// размер числа, run_siqs крутит батчи до готовности матрицы, factor<T> —
// точка входа для конвейера.
// ---------------------------------------------------------------------------

#include "common.hpp"
#include "relation.hpp"
#include "worker_pool.hpp"
#include "large_prime_graph.hpp"
#include "gauss.hpp"
#include "config.hpp"
#include "poly_batch.hpp"

namespace bignum::siqs_xor {

// --- Автоподбор параметров --------------------------------------------------

struct AutoParams
{
    u32 prime_bound;
    int block_bits;
    int num_a_factors;
};

inline AutoParams auto_tune_siqs(const Big &M)
{
    const int bits = static_cast<int>(M.bit_width());
    const double ln_M = bits * 0.6931471805599453;
    const double ln_ln_M = std::log(ln_M);

    // Граница факторной базы. Свободный член подобран по замерам: на 160/192/205/
    // 224/256 битах оптимум лежит примерно втрое ниже значения, которое давала
    // прежняя формула (+3.3). Меньше база — больше полиномов, но каждый дешевле:
    // при +2.0 время падает в 1.3…2.6 раза, а матрица Гаусса получается меньше
    // (см. таблицу производительности в README).
    const double exponent = 0.357 * std::sqrt(ln_M * ln_ln_M) + 2.0;
    const double raw_bound = std::exp(exponent);

    u32 prime_bound = static_cast<u32>(std::round(raw_bound / 10000.0) * 10000.0);
    const bool base_at_floor = prime_bound < 50000U;
    prime_bound = std::max(50000U, std::min(prime_bound, 8000000U));

    // Ширина блока. На 160–192 битах замеры дают 17 бит (B = 65536) против 16,
    // на 128 — наоборот, 16: слишком широкий интервал раздувает |Q(x)| и
    // соотношений становится меньше, чем экономит просев.
    int block_bits = 16;
    if (bits >= 150)
        block_bits = 17;
    if (bits >= 270)
        block_bits = 18;

    const double half_B = static_cast<double>(1U << (block_bits - 1));
    const double log2_target_A = 0.5 * (bits + 1.0) - std::log2(half_B);
    const double log2_avg_q = std::log2(static_cast<double>(prime_bound)) * 0.62;
    int s = static_cast<int>(std::round(log2_target_A / log2_avg_q));
    // Если база упёрлась в нижний предел, годных q заметно меньше, чем предсказывает
    // формула, и лишний множитель в A только портит полиномы: на 144–176 битах
    // s на единицу меньше даёт 1.2–1.6 раза (например, 160 бит: 0.62 → 0.40 с).
    if (base_at_floor)
        --s;
    // Потолок 10: на 300 битах формула просит 10.7, но 11 множителей дают
    // измеримо худший результат (1546 с против 1508 с) — выигрыш от вдвое более
    // редкой подготовки A не покрывает ухудшения самих полиномов.
    s = std::max(5, std::min(s, 10));

    return {prime_bound, block_bits, s};
}

inline u32 find_best_multiplier(const Big &M, const std::vector<u32> &small_primes)
{
    u32 best_k = 1;
    double best_score = -1e9;
    static const u32 cands[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21,
                                22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 41,
                                42, 43, 46, 47};
    for (u32 k : cands)
    {
        // kM должно помещаться в 512 бит без переполнения.
        if (M.bit_width() + static_cast<u32>(std::bit_width(k)) > 512)
            continue;
        const Big kM = M * Big(k);
        double score = -0.5 * std::log(static_cast<double>(k));
        const u64 mod8 = kM.low_word() & 7ULL;
        if (mod8 == 1)
            score += 2.0 * std::log(2.0);
        else if (mod8 == 5)
            score += 1.0 * std::log(2.0);
        for (size_t i = 1; i < std::min<size_t>(40, small_primes.size()); ++i)
        {
            const u32 p = small_primes[i];
            const u64 rem = kM.mod_u64(p);
            if (rem == 0)
                score += std::log(static_cast<double>(p)) / p;
            else if (pow_mod64(rem, (p - 1) / 2, p) == 1)
                score += 2.0 * std::log(static_cast<double>(p)) / (p - 1);
        }
        if (score > best_score)
        {
            best_score = score;
            best_k = k;
        }
    }
    return best_k;
}

inline std::vector<u32> sieve_primes(u32 limit)
{
    std::vector<uint8_t> is_p(static_cast<size_t>(limit) + 1, 1);
    is_p[0] = 0;
    if (limit >= 1)
        is_p[1] = 0;
    for (u32 i = 2; static_cast<u64>(i) * i <= limit; ++i)
        if (is_p[i])
            for (u64 j = static_cast<u64>(i) * i; j <= limit; j += i)
                is_p[static_cast<size_t>(j)] = 0;
    std::vector<u32> primes;
    for (u32 i = 2; i <= limit; ++i)
        if (is_p[i])
            primes.push_back(i);
    return primes;
}

// --- Основной запуск --------------------------------------------------------

inline std::optional<Big> run_siqs(const Big &M, u32 prime_bound, int block_bits,
                                   int num_a_factors,
                                   const std::function<void(const std::string &)> &status_cb,
                                   const std::function<void(double, const std::string &)> &progress_cb)
{
    const auto t0 = std::chrono::steady_clock::now();
    auto notify = [&](const std::string &msg)
    {
        if (status_cb)
            status_cb(msg);
    };
    auto report = [&](double fraction, const std::string &stage)
    {
        if (progress_cb)
            progress_cb(fraction, stage);
    };

    const std::vector<u32> raw_primes = sieve_primes(prime_bound);
    const u32 k = find_best_multiplier(M, raw_primes);
    const Big kM = M * Big(k);
    report(0.02, "подготовка решета");

    SIQSConfig cfg;
    cfg.M = M;
    cfg.kM = kM;
    cfg.B_size = 1U << block_bits;
    cfg.half_B = cfg.B_size / 2;
    cfg.s = num_a_factors;
    cfg.small_prime_skip = 17;

    for (u32 p : raw_primes)
    {
        if (p == 2)
        {
            cfg.prime_base.push_back(2);
            const u32 r = static_cast<u32>(kM.low_word() & 1ULL);
            cfg.root1.push_back(r);
            cfg.root2.push_back(r);
            cfg.log_primes.push_back(1);
        }
        else
        {
            const u64 rem = kM.mod_u64(p);
            if (rem == 0 || pow_mod64(rem, (p - 1) / 2, p) == 1)
            {
                const u32 r = tonelli_shanks(rem, p);
                cfg.prime_base.push_back(p);
                cfg.root1.push_back(r);
                cfg.root2.push_back((r != p - r) ? (p - r) : r);
                cfg.log_primes.push_back(
                    static_cast<uint8_t>(std::llround(std::log2(static_cast<double>(p)))));
            }
        }
    }

    cfg.num_primes = static_cast<u32>(cfg.prime_base.size());
    if (cfg.num_primes < 64)
        return std::nullopt;

    for (u32 i = 0; i < cfg.num_primes; ++i)
        cfg.prime_to_idx[cfg.prime_base[i]] = i;
    cfg.last_prime = cfg.prime_base.back();
    cfg.min_2lp_comp = cfg.last_prime * cfg.last_prime;
    cfg.lp_bound = static_cast<u64>(prime_bound) * 30ULL;
    cfg.lp2_bound = cfg.lp_bound * cfg.lp_bound;

    cfg.split_idx = cfg.num_primes;
    for (u32 i = 0; i < cfg.num_primes; ++i)
    {
        if (cfg.prime_base[i] >= cfg.B_size)
        {
            cfg.split_idx = i;
            break;
        }
    }

    UInt<1024> two_kM(kM);
    two_kM.shl(1);
    UInt<1024> target_A = ops::isqrt<1024>(two_kM);
    target_A.divrem_u64(cfg.half_B);
    const double log2_A = static_cast<double>(target_A.bit_width());
    const u32 ideal_q = static_cast<u32>(
        std::max(50.0, std::pow(2.0, log2_A / static_cast<double>(cfg.s))));

    u32 center_idx = cfg.num_primes / 2;
    for (u32 i = 0; i < cfg.num_primes; ++i)
    {
        if (cfg.prime_base[i] >= ideal_q)
        {
            center_idx = i;
            break;
        }
    }
    u32 pool_start = 15;
    if (cfg.num_primes > 200)
    {
        const u32 lo = (center_idx > 60) ? (center_idx - 60) : 15U;
        pool_start = std::max(15U, std::min(lo, cfg.num_primes - 130));
    }
    for (u32 i = pool_start; i < std::min(cfg.num_primes - 1, pool_start + 120); ++i)
        cfg.pool_indices.push_back(i);

    unsigned threads = cpu::worker_thread_count();
    if (threads == 0)
        threads = 1;

    notify("SIQS: " + std::to_string(M.bit_width()) + " бит, множитель k="
         + std::to_string(k) + ", база: " + std::to_string(cfg.num_primes)
         + " простых до " + std::to_string(prime_bound)
         + ", полиномов на A: " + std::to_string(1ULL << (cfg.s - 1))
         + ", потоков: " + std::to_string(threads));
    report(0.05, "сбор соотношений");

    MarkowitzGaussianGF2 gauss;
    gauss.status_cb = status_cb;
    gauss.progress_cb = [&](double f, const std::string &t)
    {
        report(0.86 + 0.12 * f, t);
    };
    LargePrimeGraph lp_graph(M);
    u64 total_polys = 0;
    u64 batch_seed = 1000;
    const int polys_per_A = static_cast<int>(1ULL << (cfg.s - 1));
    // Ориентир для полосы прогресса: для готовности матрицы нужно примерно
    // столько соотношений, сколько активных столбцов (простых базы).
    const double relation_target = static_cast<double>(std::max<u32>(cfg.num_primes, 64));

    // Нарастающая порция полиномов на поток: начинаем с двух и удваиваем, пока
    // матрица не готова. Ранняя остановка внутри батча (raw_found) не даёт
    // порции разрастись сверх нужного, но стартовый размер оставлен прежним:
    // он задаёт траекторию A-наборов, и менять её без нужды не стоит.
    int poly_limit = std::min(polys_per_A, 2);

    WorkerPool pool(threads);
    std::atomic<u64> raw_found{0};
    SieveStats total_stats;
    double t_lp_merge = 0.0;   // серийная обработка сырых соотношений
    double t_gauss_total = 0.0; // все попытки исключения
    double t_batch_wall = 0.0;  // чистое время батчей (ожидание самого медленного)

    while (true)
    {
        if (stop_requested())
            return std::nullopt;

        // Порог ранней остановки батча: сырых соотношений нужно примерно
        // втрое больше, чем строк матрицы, потому что большая часть уходит в
        // LP-циклы. Ниже порога батч не трогаем — иначе клип на средних
        // размерах сдвигает траекторию A-наборов без пользы.
        const u64 raw_stop_at = static_cast<u64>(3.0 * relation_target);
        raw_found.store(0, std::memory_order_relaxed);

        std::vector<BatchResult> batch_results(threads);
        const auto t_batch0 = std::chrono::steady_clock::now();
        pool.run([&](unsigned t)
        {
            batch_results[t] = process_one_A_batch(cfg, batch_seed + t + 1,
                                                   &raw_found, raw_stop_at, poly_limit);
        });
        t_batch_wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_batch0).count();
        batch_seed += threads;
        poly_limit = std::min(polys_per_A, poly_limit * 2);

        const auto t_lp0 = std::chrono::steady_clock::now();
        for (unsigned t = 0; t < threads; ++t)
        {
            total_stats.add(batch_results[t].stats);
            for (const auto &r : batch_results[t].relations)
            {
                FullRelation full_rel;
                if (lp_graph.add_edge_and_extract_cycle(r, full_rel))
                    gauss.add_relation(std::move(full_rel));
            }
        }
        t_lp_merge += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_lp0).count();
        total_polys = total_stats.polynomials;

        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - t0).count();
        auto sec = [](double v)
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.1f", v);
            return std::string(buf);
        };
        notify("SIQS: " + std::to_string(static_cast<long long>(elapsed)) + " с, полиномов "
             + std::to_string(total_polys) + ", соотношений: "
             + std::to_string(gauss.relations.size()) + " ("
             + std::to_string(total_stats.rel_full) + " полных / "
             + std::to_string(total_stats.rel_1lp) + " 1LP / "
             + std::to_string(total_stats.rel_2lp) + " 2LP), простых (>=2 вх.): "
             + std::to_string(gauss.multi_hit_primes) + " (синглетонов: "
             + std::to_string(gauss.singleton_primes) + ", циклов 2LP: "
             + std::to_string(lp_graph.cycles_found) + ")");
        notify("SIQS: выживших " + std::to_string(total_stats.survivors)
             + ", попаданий крупных простых " + std::to_string(total_stats.large_prime_hits)
             + " (в выживших " + std::to_string(total_stats.hits_in_survivors) + "), отброшено: "
             + std::to_string(total_stats.dropped_big_cofactor) + " по границе / "
             + std::to_string(total_stats.dropped_nonsplit) + " не разложилось; фазы: A "
             + sec(total_stats.t_prep) + " с, сито " + sec(total_stats.t_sieve)
             + " с, скан " + sec(total_stats.t_scan) + " с, разложение "
             + sec(total_stats.t_factor) + " с, LP " + sec(t_lp_merge)
             + " с, Гаусс " + sec(t_gauss_total) + " с, батчи " + sec(t_batch_wall) + " с");

        // Полоса сбора соотношений: 5%..85% по доле от ориентира. Если матрица
        // ещё не готова, полоса просто упирается в верхнюю границу — дальнейшие
        // батчи продолжают уточнять матрицу, а не «съедают» прогресс.
        report(0.05 + 0.80 * std::min(1.0, static_cast<double>(gauss.relations.size()) / relation_target),
               "сбор соотношений");

        if (!gauss.ready_to_check_rank())
            continue;

        report(0.86, "решение матрицы GF(2)");
        const auto tg0 = std::chrono::steady_clock::now();
        const auto deps = gauss.solve_pruned_matrix();
        t_gauss_total += std::chrono::duration<double>(std::chrono::steady_clock::now() - tg0).count();
        for (const auto &dep : deps)
        {
            if (stop_requested())
                return std::nullopt;

            Big LHS(1);
            std::unordered_map<u32, u32> total_exps;
            for (u32 idx_r : dep)
            {
                LHS = Big::mul_mod(LHS, gauss.relations[idx_r].X_val, M);
                for (const auto &fe : gauss.relations[idx_r].exps)
                    total_exps[fe.idx] += fe.cnt;
            }

            Big RHS(1);
            for (const auto &kv : total_exps)
            {
                const u32 p_val = cfg.prime_base[kv.first];
                const Big p_pow = pow_mod_big(Big(p_val), kv.second / 2, M);
                RHS = Big::mul_mod(RHS, p_pow, M);
            }

            const Big diff = (LHS >= RHS) ? (LHS - RHS) : (RHS - LHS);
            Big sum = LHS + RHS;
            if (sum >= M)
                sum.sub(M);

            const Big candidates[2] = {diff, sum};
            for (const Big &cand : candidates)
            {
                const Big g = ops::gcd(cand, M);
                if (g > Big(1) && g < M)
                {
                    const auto t1 = std::chrono::steady_clock::now();
                    const double total_s = std::chrono::duration<double>(t1 - t0).count();
                    const double gauss_s = t_gauss_total;
                    notify("SIQS: успех за " + std::to_string(total_s) + " с (Гаусс: "
                         + std::to_string(gauss_s) + " с, LP: "
                         + std::to_string(t_lp_merge) + " с), полиномов: "
                         + std::to_string(total_polys) + ", делитель: " + g.toString());
                    report(1.0, "делитель найден");
                    return g;
                }
            }
        }
    }
}

// Точка входа для конвейера: факторизует n, возвращает нетривиальный делитель.
// progress_cb (необязательный) получает долю 0..1 внутри решета и название
// текущей фазы — конвейер отображает её в общую полосу прогресса.
template <typename T>
std::optional<T> factor(const T &n,
                        std::function<void(const std::string &)> status_cb = nullptr,
                        std::function<void(double, const std::string &)> progress_cb = nullptr)
{
    Big M(n);
    if (M <= Big(3))
        return std::nullopt;
    if ((M.low_word() & 1ULL) == 0)
        return T(2);

    const AutoParams params = auto_tune_siqs(M);
    if (status_cb)
        status_cb("SIQS: автоподбор: prime_bound=" + std::to_string(params.prime_bound)
                  + ", block_bits=" + std::to_string(params.block_bits)
                  + ", num_a_factors=" + std::to_string(params.num_a_factors));

    const auto divisor = run_siqs(M, params.prime_bound, params.block_bits,
                                  params.num_a_factors, status_cb, progress_cb);
    if (!divisor.has_value())
        return std::nullopt;
    return T(*divisor);
}

} // namespace bignum::siqs_xor
