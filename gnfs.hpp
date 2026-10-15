#pragma once

#include "siqs.hpp"
#include "u128.hpp"
#include "u128_utils.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bignum::gnfs {

struct RationalPrime
{
    uint32_t p;
    uint8_t log_p;
};

struct AlgebraicPrime
{
    uint32_t p;
    uint32_t r; // Корень f(r) == 0 (mod p)
    uint8_t log_p;
};

struct QuadraticCharacter
{
    uint32_t q;
    uint32_t s; // Корень f(s) == 0 (mod q), f'(s) != 0 (mod q)
};

// Параметры решета, масштабируемые по разрядности числа.
// Раньше все значения были зафиксированы под 128-битные числа, из-за чего
// GNFS для чисел > ~180 бит не находил НИ ОДНОГО соотношения: типовые нормы
// вырастали до 2^140...2^180, а гладкость по базе из простых < 5000 для них
// имеет вероятность ~1e-13 на сторону (подробности в комментарии у таблицы).
struct GnfsParams
{
    uint32_t d = 4;          // степень алгебраического полинома
    size_t rfb_size = 3000;  // размер рациональной факторной базы
    size_t afb_size = 5000;  // размер алгебраической факторной базы
    size_t qcb_size = 16;    // число квадратичных характеров
    int sieve_A = 65536;     // |a| <= sieve_A
    uint32_t max_b = 60000;  // предел по b (ранее жёсткий MAX_B)
    uint64_t lp_mult = 250;  // граница large prime = p_max * lp_mult
};

// Подбор параметров по разрядности. Реальный рабочий диапазон этой реализации
// (линейное решето, base-m полином, плотный Гаусс) — до ~150 бит включительно.
// Замеры на полупростых p×q (строк ожидается / нужно в оценке достижимости):
//   135 бит — 6.2e4 / 15285 (проходит),  150 бит — 2.1e4 / 15285 (проходит),
//   152 бит — 5.7e3 / 12205 (отказ),    200 бит — 496 / 22362,
//   320 бит — 0.46 / 41658. Дефицит растёт экспоненциально, и закрыть его
// подстройкой тиров нельзя: нужны факторные базы 2^25..2^30, решетчатое решето
// (special-q) и разреженная матрица с блочным Wiedemann, т.е. другая
// реализация. Ресурсы (ядра/память) тут не помогают — не хватает выхода на пару.
// Поэтому выше 150 бит возвращаем nullopt: раньше тиры до 400 бит лишь
// удлиняли отказ — построение факторных баз с поиском корней стоило 10-54 с,
// после чего решето всё равно не запускалось.
inline std::optional<GnfsParams> select_gnfs_params(uint32_t bits)
{
    GnfsParams p;
    if (bits <= 150) {
        // Проверено практикой (128-битные полупростые). Базы 6000/9000 вместо
        // 3000/5000: больший предел гладкости заметно повышает выход на
        // 135-150 битах (u нормы падает ~ на 0.3-0.5, rho — в 2-3 раза),
        // запас по b (300000) покрывает «невезучие» числа.
        p.d = 4;  p.rfb_size = 6000;  p.afb_size = 9000;  p.qcb_size = 48;
        p.sieve_A = 65536; p.max_b = 300000; p.lp_mult = 400;
    } else {
        return std::nullopt;
    }
    return p;
}

// Верхнее простое факторной базы по её размеру (p_k-е по счёту простое).
// Нужно только для предварительной оценки достижимости ДО построения баз:
// решето Эратосфена до оценки p_k стоит миллисекунды, тогда как построение AFB
// с поиском корней по каждому простому — десятки секунд.
inline uint64_t estimate_base_max_prime(size_t entries)
{
    if (entries == 0)
        return 2;
    const double n = static_cast<double>(entries);
    // Верхняя оценка Россора—Шонфельда: p_n < n * (ln n + ln ln n) для n >= 6.
    double bound = (n >= 6.0) ? n * (std::log(n) + std::log(std::log(n))) : 16.0;
    if (bound < 16.0)
        bound = 16.0;
    auto ps = ::u128::utils::primes(static_cast<unsigned>(bound) + 1);
    const size_t idx = std::min(entries, ps.size());
    return ps[idx - 1];
}

// log10 rho(u) для функции Дикмана по таблице с линейной интерполяцией.
// Точности хватает для порядковой оценки «задача достижима / нет».
inline double log10_dickman(double u)
{
    static const double table[][2] = {
        {1.0, 0.0},    {2.0, -0.513},  {3.0, -1.313},  {4.0, -2.309},
        {5.0, -3.459}, {6.0, -4.708},  {7.0, -6.022},  {8.0, -7.391},
        {9.0, -8.806}, {10.0, -10.27}, {12.0, -13.35}, {15.0, -18.10},
        {20.0, -26.47}, {25.0, -35.43}, {30.0, -44.73}, {40.0, -64.3},
        {50.0, -84.7},
    };
    if (u <= 1.0)
        return 0.0;
    const size_t n = sizeof(table) / sizeof(table[0]);
    if (u >= table[n - 1][0])
        return -u * std::log10(u);
    for (size_t i = 1; i < n; ++i) {
        if (u <= table[i][0]) {
            double t = (u - table[i - 1][0]) / (table[i][0] - table[i - 1][0]);
            return table[i - 1][1] + t * (table[i][1] - table[i - 1][1]);
        }
    }
    return 0.0;
}

template<typename T>
struct GNFSRelation
{
    int64_t a = 0;
    uint32_t b = 0;
    bool rat_neg = false;
    bool alg_neg = false;
    std::vector<uint32_t> rat_idx;
    std::vector<uint16_t> rat_pow;
    std::vector<uint32_t> alg_idx;
    std::vector<uint16_t> alg_pow;

    uint64_t large_rat_p = 1;
    uint64_t large_alg_p = 1;
    uint32_t large_alg_r = 0;

    bool has_pair = false;
    int64_t a2 = 0;
    uint32_t b2 = 0;
};

template<typename T>
T nth_root(const T &n, uint32_t d)
{
    if (n <= T(1) || d <= 1)
        return n;
    uint32_t bits = n.bit_width();
    uint32_t root_bits = (bits + d - 1) / d;

    T low = T(1) << ((root_bits > 1) ? (root_bits - 1) : 0);
    T high = (root_bits < T::BITS) ? (T(1) << root_bits) : ~T(0);
    T ans = low;

    while (low <= high) {
        T diff = high - low;
        T mid = low + (diff >> 1);

        T pow_val = T(1);
        bool overflow = false;
        for (uint32_t i = 0; i < d; ++i) {
            if (mid > T(0) && pow_val > ::bignum::siqs::div_q(~T(0), mid)) {
                overflow = true;
                break;
            }
            pow_val = pow_val * mid;
        }

        if (!overflow && pow_val <= n) {
            ans = mid;
            low = mid + T(1);
        } else {
            if (mid == T(0))
                break;
            high = mid - T(1);
        }
    }
    return ans;
}

//

template<typename T>
uint32_t eval_poly_mod(const std::vector<T> &coeffs, uint32_t x, uint32_t p)
{
    uint64_t res = 0;
    for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i) {
        uint32_t c_mod = ::bignum::siqs::rem_u32(coeffs[i], p);
        res = (res * x + c_mod) % p;
    }
    return static_cast<uint32_t>(res);
}

template<typename T>
uint32_t eval_poly_deriv_mod(const std::vector<T> &coeffs, uint32_t x, uint32_t p)
{
    uint64_t res = 0;
    for (int i = static_cast<int>(coeffs.size()) - 1; i >= 1; --i) {
        uint64_t c_mod = ::bignum::siqs::rem_u32(coeffs[i], p);
        uint64_t term = (c_mod * static_cast<uint64_t>(i)) % p;
        res = (res * x + term) % p;
    }
    return static_cast<uint32_t>(res);
}

template<typename T>
T eval_homogeneous_poly(const std::vector<T> &coeffs, int64_t a, uint32_t b, bool &is_neg)
{
    size_t d = coeffs.size() - 1;
    T abs_a = T(static_cast<uint64_t>(std::abs(a)));
    T b_val = T(b);

    std::vector<T> a_pows(d + 1, T(1));
    std::vector<T> b_pows(d + 1, T(1));
    for (size_t i = 1; i <= d; ++i) {
        a_pows.at(i) = a_pows.at(i - 1) * abs_a;
        b_pows.at(i) = b_pows.at(i - 1) * b_val;
    }

    T pos_sum = T(0);
    T neg_sum = T(0);

    // F(a, b) = sum_{i=0}^d c_i * a^i * b^(d - i)
    for (size_t i = 0; i <= d; ++i) {
        T term = coeffs.at(i) * a_pows.at(i) * b_pows.at(d - i);
        bool term_neg = (a < 0) && ((i & 1) != 0);
        if (term_neg) {
            neg_sum = neg_sum + term;
        } else {
            pos_sum = pos_sum + term;
        }
    }

    if (pos_sum >= neg_sum) {
        is_neg = false;
        return pos_sum - neg_sum;
    } else {
        is_neg = true;
        return neg_sum - pos_sum;
    }
}

inline void merge_factor_lists(const std::vector<uint32_t> &idx1,
                               const std::vector<uint16_t> &pow1,
                               const std::vector<uint32_t> &idx2,
                               const std::vector<uint16_t> &pow2,
                               std::vector<uint32_t> &out_idx,
                               std::vector<uint16_t> &out_pow)
{
    std::vector<std::pair<uint32_t, uint16_t>> combined;
    combined.reserve(idx1.size() + idx2.size());
    for (size_t i = 0; i < idx1.size(); ++i)
        combined.emplace_back(idx1[i], pow1[i]);
    for (size_t i = 0; i < idx2.size(); ++i)
        combined.emplace_back(idx2[i], pow2[i]);

    std::sort(combined.begin(), combined.end(), [](const auto &x, const auto &y) {
        return x.first < y.first;
    });

    out_idx.clear();
    out_pow.clear();
    for (const auto &kv : combined) {
        if (!out_idx.empty() && out_idx.back() == kv.first) {
            out_pow.back() += kv.second;
        } else {
            out_idx.push_back(kv.first);
            out_pow.push_back(kv.second);
        }
    }
}

// Умножение полиномов по модулю монического f(x) и простого числа q
inline std::vector<uint32_t> poly_mul_mod_q(const std::vector<uint32_t> &A,
                                            const std::vector<uint32_t> &B,
                                            const std::vector<uint32_t> &f_mod,
                                            uint32_t q)
{
    size_t d = f_mod.size() - 1;
    std::vector<uint64_t> prod(2 * d - 1, 0);
    for (size_t i = 0; i < d; ++i) {
        if (A[i] == 0)
            continue;
        for (size_t j = 0; j < d; ++j) {
            prod[i + j] = (prod[i + j] + 1ULL * A[i] * B[j]) % q;
        }
    }
    for (int i = static_cast<int>(2 * d - 2); i >= static_cast<int>(d); --i) {
        uint64_t lead = prod[i] % q;
        if (lead == 0)
            continue;
        for (size_t j = 0; j < d; ++j) {
            uint64_t sub = (lead * f_mod[j]) % q;
            prod[i - d + j] = (prod[i - d + j] + q - sub) % q;
        }
    }
    std::vector<uint32_t> res(d);
    for (size_t i = 0; i < d; ++i)
        res[i] = static_cast<uint32_t>(prod[i] % q);
    return res;
}

inline std::vector<uint32_t> poly_pow_mod_q(std::vector<uint32_t> base,
                                            uint64_t exp,
                                            const std::vector<uint32_t> &f_mod,
                                            uint32_t q)
{
    size_t d = f_mod.size() - 1;
    std::vector<uint32_t> res(d, 0);
    res[0] = 1;
    while (exp > 0) {
        if (exp & 1)
            res = poly_mul_mod_q(res, base, f_mod, q);
        base = poly_mul_mod_q(base, base, f_mod, q);
        exp >>= 1;
    }
    return res;
}

// --- Быстрая проверка неприводимости по модулю q ---
// Полный перебор квадратных делителей стоит O(q^2) (для q ~ 65000 — миллиарды
// операций), поэтому для подбора КРУПНОГО q (больше бит на цифру базы => короче
// K => подъём Ньютона O(K^2) заметно дешевле) используем критерий Фробениуса:
// f степени d неприводимо над F_q <=> x^(q^d) == x (mod f) и
// gcd(f, x^(q^(d/ℓ)) - x) == 1 для всех простых делителей ℓ степени d.
inline void poly_trim_mod_q(std::vector<uint32_t> &v)
{
    while (!v.empty() && v.back() == 0)
        v.pop_back();
}

// НОД двух полиномов по модулю простого q (евклид; степени tiny, O(d^2))
inline std::vector<uint32_t> poly_gcd_mod_q(std::vector<uint32_t> a,
                                            std::vector<uint32_t> b,
                                            uint32_t q)
{
    poly_trim_mod_q(a);
    poly_trim_mod_q(b);
    while (!b.empty()) {
        const uint32_t inv_lead = ::bignum::siqs::mod_inv32(b.back(), q);
        std::vector<uint32_t> r = a;
        while (!r.empty() && r.size() >= b.size()) {
            const size_t shift = r.size() - b.size();
            const uint32_t coef = static_cast<uint32_t>(1ULL * r.back() * inv_lead % q);
            if (coef != 0) {
                for (size_t i = 0; i < b.size(); ++i) {
                    const uint64_t sub = 1ULL * b[i] * coef % q;
                    r[shift + i] = static_cast<uint32_t>((r[shift + i] + q - sub) % q);
                }
            }
            poly_trim_mod_q(r);
        }
        a = std::move(b);
        b = std::move(r);
    }
    return a;
}

inline bool is_irreducible_mod_q_fast(const std::vector<uint32_t> &f_mod, uint32_t q)
{
    const size_t d = f_mod.size() - 1;
    if (d < 2)
        return false;

    // Какие k = d/ℓ нужны для проверок НОД (ℓ — простой делитель d)
    std::vector<size_t> need_k;
    for (size_t pp = 2; pp <= d; ++pp) {
        if (d % pp == 0) {
            bool is_prime = true;
            for (size_t t = 2; t * t <= pp; ++t) {
                if (pp % t == 0) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime)
                need_k.push_back(d / pp);
        }
    }

    std::vector<uint32_t> g(d, 0);
    if (d >= 2)
        g[1] = 1; // x mod (f, q)
    for (size_t k = 1; k <= d; ++k) {
        g = poly_pow_mod_q(g, q, f_mod, q); // g <- g^q: g = x^(q^k) mod (f, q)
        if (k == d) {
            for (size_t i = 0; i < d; ++i) {
                const uint32_t want = (i == 1) ? 1u : 0u;
                if (g[i] != want)
                    return false; // x^(q^d) != x
            }
        } else if (std::find(need_k.begin(), need_k.end(), k) != need_k.end()) {
            std::vector<uint32_t> diff = g;
            diff[1] = (diff[1] + q - 1) % q; // g - x
            auto gd = poly_gcd_mod_q(f_mod, diff, q);
            poly_trim_mod_q(gd);
            const bool is_one = gd.empty() || (gd.size() == 1 && gd[0] != 0);
            if (!is_one)
                return false; // есть делитель степени 1..ℓ
        }
    }
    return true;
}

//
// --- Основной класс Решета Числового Поля (GNFS) ---

template<typename T>
class GNFSSolver
{
    T m_N;
    uint32_t m_d = 3;
    T m_m;
    std::vector<T> m_coeffs;

    std::vector<RationalPrime> m_RFB;
    std::vector<AlgebraicPrime> m_AFB;
    std::vector<QuadraticCharacter> m_QCB;

    int m_sieve_A = 65536;
    size_t m_relations_needed = 0;
    uint64_t m_L_rat_max = 0;
    uint64_t m_L_alg_max = 0;

    // Масштабируемые параметры (см. GnfsParams)
    uint32_t m_max_b = 60000;
    uint64_t m_lp_mult = 250;
    uint8_t m_rat_thresh = 28;
    uint8_t m_alg_thresh = 34;
    std::chrono::steady_clock::time_point m_sieve_start{};
    // Последний вывод статуса (под m_rel_mutex) — троттлинг.
    std::chrono::steady_clock::time_point m_last_status{};

    // Диагностика решета (для отчёта после прогона)
    std::atomic<uint64_t> m_stat_hits{0};        // пар прошли пороги
    std::atomic<uint64_t> m_stat_pre_rat_fail{0}; // отсеял предфильтр (рац. сторона)
    std::atomic<uint64_t> m_stat_pre_alg_fail{0}; // отсеял предфильтр (алг. сторона)
    std::atomic<uint64_t> m_stat_rat_fail{0};     // не разложилась по RFB
    std::atomic<uint64_t> m_stat_afb_fail{0};     // не разложилась по AFB
    std::atomic<uint64_t> m_stat_success{0};      // принято соотношений

    std::mutex m_rel_mutex;
    std::vector<GNFSRelation<T>> m_relations;
    std::unordered_set<uint64_t> m_seen_ab;

    // Таблицы для спаривания одиночных больших простых:
    std::unordered_map<uint64_t, GNFSRelation<T>> m_rat_partials;
    std::unordered_map<uint64_t, GNFSRelation<T>> m_alg_partials;
    // Спаривание соотношений с большими простыми С ОБЕИХ сторон:
    // ключ — (ratLP, algLP<<32 | корень). Обе пары взаимно уничтожаются.
    std::map<std::pair<uint64_t, uint64_t>, GNFSRelation<T>> m_dblp_partials;

    size_t m_full_relations_count = 0;
    size_t m_combined_count = 0;
    std::atomic<bool> m_done{false};

public:
    explicit GNFSSolver(const T &n)
        : m_N(n)
    {}

    std::optional<T> solve(std::function<void(const std::string &)> status_cb = nullptr)
    {
        uint32_t bits = m_N.bit_width();

        // Масштабируемые параметры под разрядность. Числа крупнее рабочей
        // полосы (~150 бит) получают честный отказ мгновенно: в конвейере
        // факторизации GNFS выше 150 бит и не вызывается, поэтому пауза перед
        // отказом здесь не нужна.
        auto params_opt = select_gnfs_params(bits);
        if (!params_opt.has_value()) {
            if (status_cb) {
                status_cb("GNFS: " + std::to_string(bits)
                          + "-битное число вне пределов этой реализации GNFS "
                            "(рабочий диапазон до ~150 бит; для больших чисел "
                            "нужен, например, CADO-NFS на кластере).");
            }
            return std::nullopt;
        }
        const GnfsParams &params = *params_opt;
        m_d = params.d;
        m_sieve_A = params.sieve_A;
        m_max_b = params.max_b;
        m_lp_mult = params.lp_mult;

        if (status_cb)
            status_cb("GNFS: " + std::to_string(bits) + " бит: подбор Base-m полинома (d="
                      + std::to_string(m_d) + ")...");
        if (!select_polynomial(params)) {
            return std::nullopt;
        }

        // Пороги решета выводятся из реального m: полезная зона малых b даёт
        // нормы примерно от m (см. |b*m - a| >= m - A). Для проверенной точки
        // (128 бит, m=2^32) формула воспроизводит зафиксированные 28/34.
        {
            uint32_t m_bits = m_m.bit_width();
            m_rat_thresh = static_cast<uint8_t>(std::clamp<int>(static_cast<int>(m_bits) - 4, 16, 250));
            m_alg_thresh = static_cast<uint8_t>(std::clamp<int>(static_cast<int>(m_bits) + 2, 16, 250));
        }

        if (status_cb)
            status_cb("GNFS: Полином d=" + std::to_string(m_d) + ", m=" + m_m.toString()
                      + " (" + std::to_string(m_m.bit_width()) + " бит), пороги "
                      + std::to_string(m_rat_thresh) + "/" + std::to_string(m_alg_thresh)
                      + ", область |a|<=" + std::to_string(m_sieve_A));

        // Предварительная оценка достижимости ДО построения факторных баз:
        // поиск корней для AFB стоит 5-54 с, а для оценки нужны только полином
        // и верхние простые баз. RFB — это первые rfb_size простых, так что её
        // p_max берём точно; AFB состоит из пар (p, r) с ~1.1 корня на простое,
        // поэтому для afb_size пар хватает ~0.9*afb_size простых (меньший p_max
        // делает оценку слегка строже — в консервативную сторону).
        {
            const size_t est_cols = (1 + params.rfb_size + 1 + params.afb_size)
                                    + params.qcb_size;
            const size_t est_needed = est_cols + std::max<size_t>(64, est_cols / 64);
            const uint64_t est_rfb_max = estimate_base_max_prime(params.rfb_size);
            const uint64_t est_afb_max = estimate_base_max_prime(params.afb_size * 9 / 10);
            if (!estimate_feasibility(est_rfb_max, est_afb_max, est_needed,
                                      "Предварительная оценка (до построения баз)", 0,
                                      status_cb)) {
                return std::nullopt;
            }
        }

        if (status_cb)
            status_cb("GNFS: Построение баз (RFB " + std::to_string(params.rfb_size)
                      + " / AFB " + std::to_string(params.afb_size) + ")...");
        if (!build_factor_bases(params.rfb_size, params.afb_size, params.qcb_size)) {
            return std::nullopt;
        }

        size_t total_cols = (1 + m_RFB.size() + 1 + m_AFB.size()) + m_QCB.size();
        // Строк должно быть заметно больше числа столбцов, иначе в СЛАУ не будет
        // нетривиального вектора в ядре.
        m_relations_needed = total_cols + std::max<size_t>(64, total_cols / 64);

        // Итоговая оценка достижимости по уже построенным базам (точные p_max).
        // Если бы предварительная оценка ошиблась — честный отказ здесь.
        if (!estimate_feasibility(m_RFB.back().p, m_AFB.back().p,
                                  m_relations_needed,
                                  "Итоговая оценка по построенным базам", 5000,
                                  status_cb)) {
            return std::nullopt;
        }

        m_relations.clear();
        m_seen_ab.clear();
        m_rat_partials.clear();
        m_alg_partials.clear();
        m_dblp_partials.clear();
        m_full_relations_count = 0;
        m_combined_count = 0;
        m_stat_hits.store(0);
        m_stat_pre_rat_fail.store(0);
        m_stat_pre_alg_fail.store(0);
        m_stat_rat_fail.store(0);
        m_stat_afb_fail.store(0);
        m_stat_success.store(0);
        m_done.store(false);
        m_sieve_start = std::chrono::steady_clock::now();

        if (status_cb)
            status_cb("GNFS: Запуск решета (RFB + AFB + Large Primes), b до "
                      + std::to_string(m_max_b) + "...");
        unsigned num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 4;

        std::atomic<uint32_t> b_counter{1};
        std::vector<std::thread> workers;
        workers.reserve(num_threads);

        for (unsigned t = 0; t < num_threads; ++t) {
            workers.emplace_back(
                [this, &b_counter, status_cb]() { sieve_worker(b_counter, status_cb); });
        }

        for (auto &w : workers) {
            if (w.joinable())
                w.join();
        }

        if (::u128::Globals::LoadStop()) {
            return std::nullopt;
        }

        // Диагностика решета: сколько пар прошло пороги, сколько отсеял
        // предфильтр по логам и сколько не разложилось по базам.
        if (status_cb) {
            char buf[384];
            std::snprintf(buf, sizeof(buf),
                          "GNFS: Решето: пар прошло пороги %llu, предфильтр отсеял %llu "
                          "(рац.) / %llu (алг.), не разложилось %llu / %llu, принято %llu",
                          static_cast<unsigned long long>(m_stat_hits.load()),
                          static_cast<unsigned long long>(m_stat_pre_rat_fail.load()),
                          static_cast<unsigned long long>(m_stat_pre_alg_fail.load()),
                          static_cast<unsigned long long>(m_stat_rat_fail.load()),
                          static_cast<unsigned long long>(m_stat_afb_fail.load()),
                          static_cast<unsigned long long>(m_stat_success.load()));
            status_cb(buf);
        }

        // 1. Проверка на пустоту или недостаток соотношений перед матрицей
        if (m_relations.empty() || m_relations.size() <= total_cols) {
            if (status_cb) {
                status_cb("GNFS: Остановлено! Недостаточно строк для СЛАУ: "
                          + std::to_string(m_relations.size()) + " <= " + std::to_string(total_cols)
                          + " (чистых: " + std::to_string(m_full_relations_count)
                          + ", пар: " + std::to_string(m_combined_count) + ")");
            }
            // Задержка, чтобы успеть прочитать статус на экране / в UI
            std::this_thread::sleep_for(std::chrono::milliseconds(5000));
            return std::nullopt;
        }

        if (status_cb)
            status_cb("GNFS: Решение матрицы СЛАУ над GF(2)...");

        auto result = solve_linear_algebra(status_cb);

        // 2. Если СЛАУ завершилась неудачей, выводим статус и делаем паузу
        if (!result.has_value() && status_cb) {
            status_cb("GNFS: СЛАУ не нашла нетривиальных делителей (возврат std::nullopt)");
            std::this_thread::sleep_for(std::chrono::milliseconds(5000));
        }

        return result;
    }

private:
    // Оценка качества base-m полинома: максимум битовой длины коэффициентов
    // (доминирующий размер норм).
    uint32_t poly_score(const std::vector<T> &coeffs) const
    {
        uint32_t max_bits = 0;
        for (const auto &c : coeffs) {
            uint32_t b = c.bit_width();
            if (b > max_bits)
                max_bits = b;
        }
        return max_bits;
    }

    // Подбор Base-m полинома: перебираем m в допустимом диапазоне
    // [N^(1/d), 2^(1/d) * N^(1/d)) (тогда старший коэффициент равен 1)
    // и выбираем кандидата с наименьшей нормой коэффициентов. Степень берётся
    // только из таблицы параметров: попытка d+1 для малых чисел проигрывала —
    // коэффициенты мельче, но нормы F(a,b) ~ m*(|a|+b)^d растут со степенью.
    bool select_polynomial(const GnfsParams &params)
    {
        bool found = false;
        uint32_t best_score = UINT32_MAX;
        const uint32_t d_try = params.d;

        {
            T m0 = nth_root(m_N, d_try);
            if (m0 <= T(2))
                return false;

            // Верхняя граница диапазона: (m+1)^d > N, т.е. m < 2^(1/d) * m0.
            // Оцениваем множителем 1 + 1/d (немного больше нужного — лишние
            // кандидаты отбрасывает проверка старшего коэффициента).
            constexpr int NUM_CANDIDATES = 48;
            for (int k = 0; k < NUM_CANDIDATES; ++k) {
                // m = m0 + k * step, step ~ m0 / NUM_CANDIDATES / d
                T step = m0 / T(static_cast<uint64_t>(NUM_CANDIDATES * d_try));
                if (step == T(0))
                    step = T(1);
                T cand = m0 + step * T(static_cast<uint64_t>(k));
                if (cand <= T(2))
                    continue;

                T m_pow_d = T(1);
                bool overflow = false;
                for (uint32_t i = 0; i < d_try; ++i) {
                    if (m_pow_d > ::bignum::siqs::div_q(~T(0), cand)) {
                        overflow = true;
                        break;
                    }
                    m_pow_d = m_pow_d * cand;
                }
                if (overflow || m_pow_d > m_N)
                    continue;

                auto [top_q, top_r] = ::bignum::siqs::div_rem(m_N, m_pow_d);
                (void)top_r;
                if (top_q != T(1))
                    continue; // Старший коэффициент != 1 — не подходит

                std::vector<T> coeffs;
                T temp = m_N;
                for (uint32_t i = 0; i < d_try; ++i) {
                    auto [q, r] = ::bignum::siqs::div_rem(temp, cand);
                    coeffs.push_back(r);
                    temp = q;
                }
                coeffs.push_back(temp);

                uint32_t score = poly_score(coeffs);
                if (score < best_score) {
                    best_score = score;
                    m_d = d_try;
                    m_m = cand;
                    m_coeffs = std::move(coeffs);
                    found = true;
                }
            }
        }

        return found && m_coeffs.size() == static_cast<size_t>(m_d) + 1;
    }

    bool build_factor_bases(size_t rfb_target, size_t afb_target, size_t qcb_target)
    {
        size_t max_prime_search = std::max<size_t>(200000, afb_target * 15);
        auto all_primes = ::u128::utils::primes(max_prime_search);

        m_RFB.clear();
        m_AFB.clear();
        m_QCB.clear();

        size_t rfb_limit = std::min(rfb_target, all_primes.size());
        m_RFB.reserve(rfb_limit);
        for (size_t i = 0; i < rfb_limit; ++i) {
            uint32_t p = all_primes[i];
            uint8_t lp = static_cast<uint8_t>(std::max(1.0, std::round(std::log2(p))));
            m_RFB.push_back({p, lp});
        }

        size_t prime_scan_limit = std::min(all_primes.size(), afb_target * 3 / 2 + 500);
        unsigned num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 4;

        std::vector<std::vector<AlgebraicPrime>> local_afb(num_threads);
        std::atomic<size_t> next_idx{0};
        std::atomic<size_t> total_found{0};

        std::vector<std::thread> workers;
        workers.reserve(num_threads);

        for (unsigned t = 0; t < num_threads; ++t) {
            workers.emplace_back([&, t]() {
                const size_t deg = m_coeffs.size() - 1;
                std::vector<uint32_t> c_mod(deg + 1);

                while (!::u128::Globals::LoadStop()) {
                    if (total_found.load(std::memory_order_relaxed) >= afb_target + 256)
                        break;

                    size_t idx = next_idx.fetch_add(32, std::memory_order_relaxed);
                    if (idx >= prime_scan_limit)
                        break;

                    size_t end_idx = std::min(idx + 32, prime_scan_limit);
                    for (size_t k = idx; k < end_idx; ++k) {
                        uint32_t p = all_primes[k];
                        uint8_t lp = static_cast<uint8_t>(std::max(1.0, std::round(std::log2(p))));

                        for (size_t i = 0; i <= deg; ++i) {
                            c_mod[i] = ::bignum::siqs::rem_u32(m_coeffs[i], p);
                        }

                        for (uint32_t r = 0; r < p; ++r) {
                            uint64_t val = c_mod[deg];
                            for (int i = static_cast<int>(deg) - 1; i >= 0; --i) {
                                val = (val * r + c_mod[i]) % p;
                            }
                            if (val == 0) {
                                local_afb[t].push_back({p, r, lp});
                                total_found.fetch_add(1, std::memory_order_relaxed);
                            }
                        }
                    }
                }
            });
        }

        for (auto &w : workers) {
            if (w.joinable())
                w.join();
        }

        if (::u128::Globals::LoadStop())
            return false;

        for (auto &vec : local_afb) {
            m_AFB.insert(m_AFB.end(), vec.begin(), vec.end());
        }
        std::sort(m_AFB.begin(), m_AFB.end(), [](const AlgebraicPrime &a, const AlgebraicPrime &b) {
            return (a.p != b.p) ? (a.p < b.p) : (a.r < b.r);
        });

        if (m_AFB.size() > afb_target) {
            m_AFB.resize(afb_target);
        }

        if (m_RFB.empty() || m_AFB.empty())
            return false;

        uint32_t max_used_p = std::max(m_RFB.back().p, m_AFB.back().p);
        auto it = std::upper_bound(all_primes.begin(), all_primes.end(), max_used_p);
        size_t q_idx = std::distance(all_primes.begin(), it);

        const size_t deg = m_coeffs.size() - 1;
        std::vector<uint32_t> c_mod(deg + 1);

        for (; q_idx < all_primes.size() && m_QCB.size() < qcb_target; ++q_idx) {
            if (::u128::Globals::LoadStop())
                return false;
            uint32_t q = all_primes[q_idx];

            for (size_t i = 0; i <= deg; ++i) {
                c_mod[i] = ::bignum::siqs::rem_u32(m_coeffs[i], q);
            }

            for (uint32_t s = 0; s < q; ++s) {
                uint64_t val = c_mod[deg];
                for (int i = static_cast<int>(deg) - 1; i >= 0; --i) {
                    val = (val * s + c_mod[i]) % q;
                }
                if (val == 0 && eval_poly_deriv_mod(m_coeffs, s, q) != 0) {
                    m_QCB.push_back({q, s});
                    if (m_QCB.size() >= qcb_target)
                        break;
                }
            }
        }

        uint64_t max_rfb_p = m_RFB.back().p;
        uint64_t max_afb_p = m_AFB.back().p;

        // Граница large prime из параметров (p_max * lp_mult); пул не раздуваем,
        // чтобы одиночные LP постоянно сталкивались в пары:
        m_L_rat_max = max_rfb_p * m_lp_mult;
        m_L_alg_max = max_afb_p * m_lp_mult;

        return true;
    }

    // Оценка достижимости: интеграл плотности гладких пар по всей области
    // решета (|a| <= A, b <= max_b) квадратурой по лог-сеткам. Равномерная
    // или «малопарная» выборка не годится: вклад дают только пары с малыми
    // |a| и b (нормы растут как m*(|a|+b)^d), а их доля мала. Считаем
    // rho(u_rat)*rho(u_alg) по узлам сетки и домножаем на ширины интервалов.
    // Если ожидание ниже нужного числа строк — честный отказ сразу, без
    // многочасового пустого прогона.
    // Верхние простые баз и нужное число строк передаются снаружи: это позволяет
    // запускать оценку ДО построения факторных баз (p_max берутся из оценки по
    // размеру базы) и экономить 10-54 с поиска корней на заведомо нереалистичных
    // числах. tag различает предварительную и итоговую оценки, refusal_delay_ms —
    // пауза перед отказом, чтобы пользователь успел прочитать сообщение (0 в
    // предварительной оценке: там важна скорость отказа).
    bool estimate_feasibility(uint64_t p_max_rfb, uint64_t p_max_afb,
                              size_t relations_needed, const char *tag,
                              int refusal_delay_ms,
                              const std::function<void(const std::string &)> &status_cb) const
    {
        const double log2_pr = std::log2(static_cast<double>(p_max_rfb));
        const double log2_pa = std::log2(static_cast<double>(p_max_afb));

        // Лог-сетки: b и |a|. Точки концентрируются у малых значений, где
        // интегранд ненулевой.
        constexpr int NUM_B = 40;
        constexpr int NUM_A = 40;
        auto make_grid = [](double hi, int n) {
            std::vector<double> g;
            double ratio = std::pow(hi, 1.0 / (n - 1));
            double v = 1.0;
            for (int i = 0; i < n; ++i) {
                g.push_back(v);
                v *= ratio;
            }
            g.back() = hi;
            return g;
        };
        std::vector<double> b_grid = make_grid(static_cast<double>(m_max_b), NUM_B);
        std::vector<double> a_grid = make_grid(static_cast<double>(m_sieve_A), NUM_A);

        double expected = 0.0;
        // Запас от больших простых (одиночный LP заметно повышает вероятность;
        // коэффициент консервативный).
        const double lp_boost = 2.0;

        for (int i = 0; i < NUM_B; ++i) {
            double b_lo = b_grid[static_cast<size_t>(i)];
            double b_hi = (i + 1 < NUM_B) ? b_grid[static_cast<size_t>(i + 1)]
                                          : static_cast<double>(m_max_b) * 1.5;
            double b_mid = 0.5 * (b_lo + b_hi);
            double width_b = b_hi - b_lo;
            uint32_t b_val = static_cast<uint32_t>(std::max(1.0, std::round(b_mid)));

            for (int j = 0; j < NUM_A; ++j) {
                double a_lo = a_grid[static_cast<size_t>(j)];
                double a_hi = (j + 1 < NUM_A) ? a_grid[static_cast<size_t>(j + 1)]
                                              : static_cast<double>(m_sieve_A) * 1.5;
                double a_mid = 0.5 * (a_lo + a_hi);
                // Оба знака a.
                double width_a = 2.0 * (a_hi - a_lo);
                int64_t a_val = std::max<int64_t>(1, static_cast<int64_t>(std::round(a_mid)));

                T bm = T(b_val) * m_m;
                T abs_a(static_cast<uint64_t>(a_val));
                T rat_norm = (abs_a >= bm) ? (abs_a - bm) : (bm - abs_a);
                if (rat_norm == T(0))
                    continue;
                bool neg = false;
                T alg_norm = eval_homogeneous_poly(m_coeffs, a_val, b_val, neg);
                if (alg_norm == T(0))
                    continue;

                double ur = static_cast<double>(rat_norm.bit_width()) / log2_pr;
                double ua = static_cast<double>(alg_norm.bit_width()) / log2_pa;
                double p_pair = std::pow(10.0, log10_dickman(ur) + log10_dickman(ua));
                expected += p_pair * width_b * width_a;
            }
        }
        expected *= lp_boost;

        // Эмпирическая калибровка: rho-оценка систематически завышает выход
        // в ~2 раза (корреляция гладкости двух норм, кратности, пороги решета).
        // Замеры: 128 бит — предсказание 4.8e4 при фактическом выходе ~2.4e4;
        // 144 бит — 1.4e4 против ~7e3.
        expected *= 0.5;

        if (status_cb) {
            char buf[320];
            std::snprintf(buf, sizeof(buf),
                          "GNFS: %s: ожидается ~%.3g строк "
                          "(нужно %zu, базы %llu / %llu)",
                          tag, expected, relations_needed,
                          static_cast<unsigned long long>(p_max_rfb),
                          static_cast<unsigned long long>(p_max_afb));
            status_cb(buf);
        }

        if (expected < static_cast<double>(relations_needed)) {
            if (status_cb) {
                char buf[448];
                std::snprintf(buf, sizeof(buf),
                              "GNFS: Задача оценена как нереалистичная для этой "
                              "реализации (ожидается ~%.3g строк, нужно %zu при "
                              "базах %llu / %llu). Число слишком велико: нужен, "
                              "например, CADO-NFS на кластере.",
                              expected, relations_needed,
                              static_cast<unsigned long long>(p_max_rfb),
                              static_cast<unsigned long long>(p_max_afb));
                status_cb(buf);
                if (refusal_delay_ms > 0)
                    std::this_thread::sleep_for(std::chrono::milliseconds(refusal_delay_ms));
            }
            return false;
        }
        return true;
    }

    // Объединение двух соотношений, разделяющих одни и те же большие простые,
    // в одну строку (LP-пары становятся полными квадратами). Поля больших
    // простых выставляет вызывающий код — они различаются по типу спаривания.
    static GNFSRelation<T> combine_pair(const GNFSRelation<T> &first,
                                        const GNFSRelation<T> &second)
    {
        GNFSRelation<T> combined;
        combined.a = first.a;
        combined.b = first.b;
        combined.has_pair = true;
        combined.a2 = second.a;
        combined.b2 = second.b;
        combined.rat_neg = (first.rat_neg != second.rat_neg);
        combined.alg_neg = (first.alg_neg != second.alg_neg);

        merge_factor_lists(first.rat_idx, first.rat_pow, second.rat_idx, second.rat_pow,
                           combined.rat_idx, combined.rat_pow);
        merge_factor_lists(first.alg_idx, first.alg_pow, second.alg_idx, second.alg_pow,
                           combined.alg_idx, combined.alg_pow);
        return combined;
    }

    void sieve_worker(std::atomic<uint32_t> &b_counter,
                      std::function<void(const std::string &)> status_cb)
    {
        const int M = m_sieve_A;
        // log2-суммы гладких кандидатов у крупных баз превышают 255,
        // поэтому храним в uint16, а не в uint8.
        std::vector<uint16_t> rat_sieve(2 * M, 0);
        std::vector<uint16_t> alg_sieve(2 * M, 0);

        // Запас предфильтра "лог решета + запас < log2 нормы": граница LP
        // плюс подстраховка на кратности (решето добавляет log p один раз).
        auto lp_bits_of = [](uint64_t v) -> uint32_t {
            return (v > 1) ? static_cast<uint32_t>(64 - std::countl_zero(v - 1)) : 0u;
        };
        const uint32_t lp_slack_rat = lp_bits_of(m_L_rat_max) + 20;
        const uint32_t lp_slack_alg = lp_bits_of(m_L_alg_max) + 20;

        const uint32_t max_b = m_max_b;
        while (!m_done.load(std::memory_order_relaxed) && !::u128::Globals::LoadStop()) {
            uint32_t b = b_counter.fetch_add(1);
            if (b > max_b)
                break;

            // Статус — не чаще раза в 700 мс на все потоки (m_rel_mutex
            // сериализует вывод; иначе статусы залипают тысячами строк).
            if (status_cb && (b % 10 == 0)) {
                std::lock_guard<std::mutex> lock(m_rel_mutex);
                auto now = std::chrono::steady_clock::now();
                if (!m_done.load(std::memory_order_relaxed)
                    && now - m_last_status > std::chrono::milliseconds(700)) {
                    m_last_status = now;
                    double secs = std::chrono::duration<double>(now - m_sieve_start).count();
                    size_t ready = m_full_relations_count + m_combined_count;
                    double rate = (secs > 1.0) ? static_cast<double>(ready) / secs : 0.0;
                    double rest = (m_relations_needed > ready)
                                      ? static_cast<double>(m_relations_needed - ready)
                                      : 0.0;
                    double eta_min = (rate > 0.0) ? (rest / rate) / 60.0 : 0.0;
                    status_cb("GNFS: Сбор соотношений: " + std::to_string(ready) + " / "
                              + std::to_string(m_relations_needed) + " (чистых: "
                              + std::to_string(m_full_relations_count) + ", пар: "
                              + std::to_string(m_combined_count) + ", b=" + std::to_string(b)
                              + ", " + std::to_string(static_cast<int>(rate)) + " строк/с, ~"
                              + std::to_string(static_cast<int>(eta_min)) + " мин)");
                }
            }

            std::fill_n(rat_sieve.begin(), 2 * M, 0);
            std::fill_n(alg_sieve.begin(), 2 * M, 0);

            // 1. Рациональное решето
            for (size_t i = 0; i < m_RFB.size(); ++i) {
                uint32_t p = m_RFB.at(i).p;
                if (p <= 2)
                    continue;
                uint32_t m_mod = ::bignum::siqs::rem_u32(m_m, p);
                uint32_t root = static_cast<uint32_t>((1ULL * (b % p) * m_mod) % p);

                int64_t rem_val = (static_cast<int64_t>(root) + M) % static_cast<int64_t>(p);
                if (rem_val < 0)
                    rem_val += p;
                int start = static_cast<int>(rem_val);

                uint8_t lp = m_RFB.at(i).log_p;
                for (int j = start; j < 2 * M; j += p) {
                    rat_sieve.at(j) += lp;
                }
            }

            // 2. Алгебраическое решето
            for (size_t i = 0; i < m_AFB.size(); ++i) {
                uint32_t p = m_AFB.at(i).p;
                if (p <= 2)
                    continue;
                uint32_t b_mod = b % p;
                if (b_mod == 0)
                    continue;

                uint32_t root = static_cast<uint32_t>((1ULL * b_mod * m_AFB.at(i).r) % p);

                int64_t rem_val = (static_cast<int64_t>(root) + M) % static_cast<int64_t>(p);
                if (rem_val < 0)
                    rem_val += p;
                int start = static_cast<int>(rem_val);

                uint8_t lp = m_AFB.at(i).log_p;
                for (int j = start; j < 2 * M; j += p) {
                    alg_sieve.at(j) += lp;
                }
            }

            constexpr int NUM_BLOCKS = 16;
            const int block_len = (2 * M) / NUM_BLOCKS;

            // Пороги выводятся из реального m (см. solve): для проверенной
            // точки 128 бит (m=2^32) получается ровно 28/34.
            const uint8_t rat_thresh = m_rat_thresh;
            const uint8_t alg_thresh = m_alg_thresh;

            for (int blk = 0; blk < NUM_BLOCKS; ++blk) {
                int j_start = blk * block_len;
                int j_end = (blk == NUM_BLOCKS - 1) ? (2 * M) : (j_start + block_len);

                for (int j = j_start; j < j_end; ++j) {
                    if (rat_sieve.at(j) >= rat_thresh && alg_sieve.at(j) >= alg_thresh) {
                        m_stat_hits++;
                        int64_t a = static_cast<int64_t>(j) - M;
                        if (a == 0)
                            continue;
                        if (((std::abs(a) | b) & 1) == 0)
                            continue;
                        if (std::gcd(std::abs(a), static_cast<int64_t>(b)) != 1)
                            continue;

                        T bm = T(b) * m_m;
                        T rat_norm;
                        bool rat_neg = false;
                        if (a >= 0) {
                            T pos_a = T(static_cast<uint64_t>(a));
                            if (pos_a >= bm) {
                                rat_norm = pos_a - bm;
                            } else {
                                rat_norm = bm - pos_a;
                                rat_neg = true;
                            }
                        } else {
                            rat_norm = bm + T(static_cast<uint64_t>(-a));
                            rat_neg = true;
                        }
                        if (rat_norm == T(0))
                            continue;

                        bool alg_neg = false;
                        T alg_norm = eval_homogeneous_poly(m_coeffs, a, b, alg_neg);
                        if (alg_norm == T(0))
                            continue;

                        // Предфильтр по логам решета: решето не видит большие
                        // простые и кратности, поэтому его сумма занижает log2
                        // нормы. Если "лог решета + граница LP + запас" всё же
                        // меньше реальной длины нормы — в норме есть крупный
                        // простой вне базы, деление бессмысленно. Отсекает
                        // подавляющее большинство кандидатов за пару сравнений.
                        if (static_cast<uint32_t>(rat_sieve.at(j)) + lp_slack_rat
                            < rat_norm.bit_width()) {
                            m_stat_pre_rat_fail++;
                            continue;
                        }
                        if (static_cast<uint32_t>(alg_sieve.at(j)) + lp_slack_alg
                            < alg_norm.bit_width()) {
                            m_stat_pre_alg_fail++;
                            continue;
                        }

                        std::vector<uint32_t> r_idx;
                        std::vector<uint16_t> r_pow;
                        uint64_t large_rat = 1;
                        if (!factor_over_rfb(rat_norm, r_idx, r_pow, large_rat)) {
                            m_stat_rat_fail++;
                            continue;
                        }

                        std::vector<uint32_t> a_idx;
                        std::vector<uint16_t> a_pow;
                        uint64_t large_alg = 1;
                        uint32_t large_alg_r = 0;
                        if (!factor_over_afb(alg_norm, a, b, a_idx, a_pow, large_alg, large_alg_r)) {
                            m_stat_afb_fail++;
                            continue;
                        }

                        uint64_t ab_key = (static_cast<uint64_t>(static_cast<uint32_t>(a + M)) << 32)
                                          | b;

                        GNFSRelation<T> rel;
                        rel.a = a;
                        rel.b = b;
                        rel.rat_neg = rat_neg;
                        rel.alg_neg = alg_neg;
                        rel.rat_idx = std::move(r_idx);
                        rel.rat_pow = std::move(r_pow);
                        rel.alg_idx = std::move(a_idx);
                        rel.alg_pow = std::move(a_pow);
                        rel.large_rat_p = large_rat;
                        rel.large_alg_p = large_alg;
                        rel.large_alg_r = large_alg_r;

                        {
                            std::lock_guard<std::mutex> lock(m_rel_mutex);
                            if (!m_seen_ab.insert(ab_key).second)
                                continue;

                            // 0. Большие простые с ОБЕИХ сторон: спариваем по
                            // точному совпадению пары (ratLP, algLP + корень) —
                            // обе составляющие взаимно уничтожаются (квадраты).
                            if (large_rat > 1 && large_alg > 1) {
                                auto key = std::make_pair(
                                    large_rat,
                                    (static_cast<uint64_t>(large_alg) << 32) | large_alg_r);
                                auto it = m_dblp_partials.find(key);
                                if (it == m_dblp_partials.end()) {
                                    m_dblp_partials.emplace(key, std::move(rel));
                                } else {
                                    GNFSRelation<T> combined = combine_pair(it->second, rel);
                                    combined.large_rat_p = large_rat; // q^2 -> +2 в X
                                    combined.large_alg_p = 1;
                                    m_relations.push_back(std::move(combined));
                                    m_dblp_partials.erase(it);
                                    m_combined_count++;
                                    m_stat_success++;
                                }
                            }
                            // 1. Полностью гладкое соотношение (0-LP)
                            else if (large_rat == 1 && large_alg == 1) {
                                m_relations.push_back(std::move(rel));
                                m_full_relations_count++;
                                m_stat_success++;
                            }
                            // 2. Рациональный Large Prime (1-LP)
                            else if (large_rat > 1) {
                                auto it = m_rat_partials.find(large_rat);
                                if (it == m_rat_partials.end()) {
                                    m_rat_partials.emplace(large_rat, std::move(rel));
                                } else {
                                    // Найдена пара! Объединяем их в одну строку:
                                    GNFSRelation<T> combined = combine_pair(it->second, rel);
                                    combined.large_rat_p
                                        = large_rat; // Квадрат L^2 взаимно уничтожится mod 2
                                    combined.large_alg_p = 1;
                                    m_relations.push_back(std::move(combined));
                                    m_rat_partials.erase(it);
                                    m_combined_count++;
                                    m_stat_success++;
                                }
                            }
                            // 3. Алгебраический Large Prime (1-LP)
                            else if (large_alg > 1) {
                                uint64_t alg_key = (large_alg << 32) | large_alg_r;
                                auto it = m_alg_partials.find(alg_key);
                                if (it == m_alg_partials.end()) {
                                    m_alg_partials.emplace(alg_key, std::move(rel));
                                } else {
                                    GNFSRelation<T> combined = combine_pair(it->second, rel);
                                    combined.large_rat_p = 1;
                                    combined.large_alg_p = large_alg;
                                    combined.large_alg_r = large_alg_r;
                                    m_relations.push_back(std::move(combined));
                                    m_alg_partials.erase(it);
                                    m_combined_count++;
                                }
                            }

                            size_t total_ready = m_full_relations_count + m_combined_count;
                            // Останавливаемся ТОЛЬКО когда готовых (гарантированно
                            // парных) строк достаточно!
                            if (total_ready >= m_relations_needed) {
                                m_done.store(true, std::memory_order_relaxed);
                                return;
                            }
                        }
                    }
                }
            }
        }
        if (b_counter.load() > m_max_b && !m_done.load()) {
            if (status_cb) {
                status_cb("GNFS: Исчерпан лимит b ("
                          + std::to_string(m_max_b) + "). Собрано строк: "
                          + std::to_string(m_relations.size()) + " / "
                          + std::to_string(m_relations_needed));
                std::this_thread::sleep_for(std::chrono::milliseconds(5000));
            }
        }
    }

    bool factor_over_rfb(T rem,
                         std::vector<uint32_t> &idx_out,
                         std::vector<uint16_t> &pow_out,
                         uint64_t &large_prime_out)
    {
        large_prime_out = 1;
        if ((::u128::utils::to_u32_low(rem) & 1) == 0) {
            uint16_t c2 = 0;
            while ((::u128::utils::to_u32_low(rem) & 1) == 0 && rem > T(0)) {
                rem >>= 1;
                c2++;
            }
            idx_out.push_back(0);
            pow_out.push_back(c2);
        }

        const T max_u64_T(std::numeric_limits<uint64_t>::max());
        bool fits_u64 = (rem <= max_u64_T);
        uint64_t rem64 = fits_u64 ? ::u128::utils::to_u64_low(rem) : 0;

        for (size_t i = 1; i < m_RFB.size(); ++i) {
            uint32_t p = m_RFB.at(i).p;
            if (fits_u64) {
                if (rem64 == 1)
                    break;
                if (1ULL * p * p > rem64) {
                    // Остаток rem64 гарантированно простой! Проверяем, есть ли он в RFB:
                    if (rem64 <= m_RFB.back().p) {
                        auto it = std::lower_bound(m_RFB.begin() + i,
                                                   m_RFB.end(),
                                                   static_cast<uint32_t>(rem64),
                                                   [](const RationalPrime &rp, uint32_t val) {
                                                       return rp.p < val;
                                                   });
                        if (it != m_RFB.end() && it->p == static_cast<uint32_t>(rem64)) {
                            idx_out.push_back(
                                static_cast<uint32_t>(std::distance(m_RFB.begin(), it)));
                            pow_out.push_back(1);
                            rem64 = 1;
                        }
                    }
                    break;
                }
                if (rem64 % p == 0) {
                    uint16_t cnt = 0;
                    do {
                        cnt++;
                        rem64 /= p;
                    } while (rem64 % p == 0);
                    idx_out.push_back(static_cast<uint32_t>(i));
                    pow_out.push_back(cnt);
                }
            } else {
                if (::bignum::siqs::rem_u32(rem, p) == 0) {
                    uint16_t cnt = 0;
                    do {
                        cnt++;
                        rem = ::bignum::siqs::div_by_u32(rem, p);
                    } while (rem > T(1) && ::bignum::siqs::rem_u32(rem, p) == 0);
                    idx_out.push_back(static_cast<uint32_t>(i));
                    pow_out.push_back(cnt);
                    if (rem <= max_u64_T) {
                        fits_u64 = true;
                        rem64 = ::u128::utils::to_u64_low(rem);
                    }
                }
            }
        }

        if (!fits_u64)
            return false;
        if (rem64 == 1)
            return true;
        if (rem64 <= m_L_rat_max && rem64 > m_RFB.back().p) {
            large_prime_out = rem64;
            return true;
        }
        return false;
    }

    bool factor_over_afb(T rem,
                         int64_t a,
                         uint32_t b,
                         std::vector<uint32_t> &idx_out,
                         std::vector<uint16_t> &pow_out,
                         uint64_t &large_prime_out,
                         uint32_t &large_r_out)
    {
        large_prime_out = 1;
        large_r_out = 0;

        const T max_u64_T(std::numeric_limits<uint64_t>::max());
        bool fits_u64 = (rem <= max_u64_T);
        uint64_t rem64 = fits_u64 ? ::u128::utils::to_u64_low(rem) : 0;

        for (size_t i = 0; i < m_AFB.size(); ++i) {
            uint32_t p = m_AFB.at(i).p;
            uint32_t r = m_AFB.at(i).r;

            if (fits_u64) {
                if (rem64 == 1)
                    break;
                if (1ULL * p * p > rem64) {
                    // Быстрый выход: если простой остаток попадает в AFB, находим его бинарным поиском
                    if (rem64 <= m_AFB.back().p) {
                        uint32_t p_rem = static_cast<uint32_t>(rem64);
                        int64_t a_mod = a % static_cast<int64_t>(p_rem);
                        if (a_mod < 0)
                            a_mod += p_rem;
                        auto it = std::lower_bound(m_AFB.begin() + i,
                                                   m_AFB.end(),
                                                   p_rem,
                                                   [](const AlgebraicPrime &ap, uint32_t val) {
                                                       return ap.p < val;
                                                   });
                        while (it != m_AFB.end() && it->p == p_rem) {
                            if ((1ULL * b * it->r) % p_rem == static_cast<uint64_t>(a_mod)) {
                                idx_out.push_back(
                                    static_cast<uint32_t>(std::distance(m_AFB.begin(), it)));
                                pow_out.push_back(1);
                                rem64 = 1;
                                break;
                            }
                            ++it;
                        }
                    }
                    break;
                }

                int64_t a_mod = a % static_cast<int64_t>(p);
                if (a_mod < 0)
                    a_mod += p;
                if ((1ULL * b * r) % p == static_cast<uint64_t>(a_mod) && (rem64 % p == 0)) {
                    uint16_t cnt = 0;
                    do {
                        cnt++;
                        rem64 /= p;
                    } while (rem64 % p == 0);
                    idx_out.push_back(static_cast<uint32_t>(i));
                    pow_out.push_back(cnt);
                }
            } else {
                int64_t a_mod = a % static_cast<int64_t>(p);
                if (a_mod < 0)
                    a_mod += p;
                if ((1ULL * b * r) % p == static_cast<uint64_t>(a_mod)
                    && ::bignum::siqs::rem_u32(rem, p) == 0) {
                    uint16_t cnt = 0;
                    do {
                        cnt++;
                        rem = ::bignum::siqs::div_by_u32(rem, p);
                    } while (rem > T(1) && ::bignum::siqs::rem_u32(rem, p) == 0);
                    idx_out.push_back(static_cast<uint32_t>(i));
                    pow_out.push_back(cnt);
                    if (rem <= max_u64_T) {
                        fits_u64 = true;
                        rem64 = ::u128::utils::to_u64_low(rem);
                    }
                }
            }
        }

        if (!fits_u64)
            return false;
        if (rem64 == 1)
            return true;
        if (rem64 <= m_L_alg_max && rem64 > m_AFB.back().p && rem64 <= 0x3FFFFFFFULL) {
            uint32_t p32 = static_cast<uint32_t>(rem64);
            uint32_t b_mod = b % p32;
            if (b_mod != 0 && std::gcd(b_mod, p32) == 1) {
                int64_t pos_a = a % static_cast<int64_t>(p32);
                if (pos_a < 0)
                    pos_a += p32;
                uint32_t b_inv = ::bignum::siqs::mod_inv32(b_mod, p32);
                large_r_out = static_cast<uint32_t>((static_cast<uint64_t>(pos_a) * b_inv) % p32);
                large_prime_out = rem64;
                return true;
            }
        }
        return false;
    }

    // Проверка неприводимости монического полинома f(x) степени d по модулю q
    bool is_irreducible_mod_q(const std::vector<uint32_t> &f_mod, uint32_t q) const
    {
        size_t d = f_mod.size() - 1;
        // 1. Проверка отсутствия корней степени 1 (для любого d)
        for (uint32_t r = 0; r < q; ++r) {
            uint64_t val = 1;
            for (int i = static_cast<int>(d) - 1; i >= 0; --i) {
                val = (val * r + f_mod[i]) % q;
            }
            if (val == 0)
                return false;
        }
        if (d <= 3)
            return true;

        // 2. Для d = 4: f(x) не должно делиться на x^2 + a1*x + a0
        if (d == 4) {
            for (uint32_t a0 = 1; a0 < q; ++a0) {
                for (uint32_t a1 = 0; a1 < q; ++a1) {
                    // Деление f(x) на (x^2 + a1*x + a0) в F_q
                    // f = x^4 + c3*x^3 + c2*x^2 + c1*x + c0
                    int64_t q1 = (static_cast<int64_t>(f_mod[3]) - a1) % q;
                    if (q1 < 0)
                        q1 += q;

                    int64_t q0 = (static_cast<int64_t>(f_mod[2]) - a0 - q1 * a1) % q;
                    if (q0 < 0)
                        q0 += q;

                    int64_t r1 = (static_cast<int64_t>(f_mod[1]) - q1 * a0 - q0 * a1) % q;
                    if (r1 < 0)
                        r1 += q;

                    int64_t r0 = (static_cast<int64_t>(f_mod[0]) - q0 * a0) % q;
                    if (r0 < 0)
                        r0 += q;

                    if (r1 == 0 && r0 == 0)
                        return false;
                }
            }
            return true;
        }

        // 3. Для d = 5 (существующая проверка)
        for (uint32_t a1 = 0; a1 < q; ++a1) {
            for (uint32_t a0 = 1; a0 < q; ++a0) {
                int64_t r4 = (static_cast<int64_t>(f_mod.at(4)) - a1) % q;
                if (r4 < 0)
                    r4 += q;
                int64_t r3 = (static_cast<int64_t>(f_mod.at(3)) - a0 - r4 * a1) % q;
                if (r3 < 0)
                    r3 += q;
                int64_t r2 = (static_cast<int64_t>(f_mod.at(2)) - r4 * a0 - r3 * a1) % q;
                if (r2 < 0)
                    r2 += q;
                int64_t r1 = (static_cast<int64_t>(f_mod.at(1)) - r3 * a0 - r2 * a1) % q;
                if (r1 < 0)
                    r1 += q;
                int64_t r0 = (static_cast<int64_t>(f_mod.at(0)) - r2 * a0) % q;
                if (r0 < 0)
                    r0 += q;
                if (r1 == 0 && r0 == 0)
                    return false;
            }
        }
        return true;
    }

    // Арифметика в системе счисления по основанию q (по модулю q^K):
    // Число хранится как std::vector<uint32_t> длины K (младшие разряды в начале).
    using RadixQ = std::vector<uint32_t>;

    static RadixQ rq_add(const RadixQ &a, const RadixQ &b, uint32_t q, size_t K)
    {
        RadixQ res(K, 0);
        uint32_t carry = 0;
        for (size_t i = 0; i < K; ++i) {
            uint32_t av = (i < a.size()) ? a.at(i) : 0;
            uint32_t bv = (i < b.size()) ? b.at(i) : 0;
            uint32_t s = av + bv + carry;
            if (s >= q) {
                res.at(i) = s - q;
                carry = 1;
            } else {
                res.at(i) = s;
                carry = 0;
            }
        }
        return res;
    }

    static RadixQ rq_sub(const RadixQ &a, const RadixQ &b, uint32_t q, size_t K)
    {
        RadixQ res(K, 0);
        int32_t borrow = 0;
        for (size_t i = 0; i < K; ++i) {
            int32_t av = (i < a.size()) ? static_cast<int32_t>(a.at(i)) : 0;
            int32_t bv = (i < b.size()) ? static_cast<int32_t>(b.at(i)) : 0;
            int32_t d = av - bv - borrow;
            if (d < 0) {
                res.at(i) = static_cast<uint32_t>(d + static_cast<int32_t>(q));
                borrow = 1;
            } else {
                res.at(i) = static_cast<uint32_t>(d);
                borrow = 0;
            }
        }
        return res;
    }

    static RadixQ rq_mul(const RadixQ &a, const RadixQ &b, uint32_t q, size_t K)
    {
        std::vector<uint64_t> acc(K, 0);
        size_t a_len = std::min(a.size(), K);
        size_t b_len = std::min(b.size(), K);
        for (size_t i = 0; i < a_len; ++i) {
            uint64_t ai = a.at(i);
            if (ai == 0)
                continue;
            size_t max_j = std::min(b_len, K - i);
            for (size_t j = 0; j < max_j; ++j) {
                acc.at(i + j) += ai * b.at(j);
            }
        }
        RadixQ res(K, 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < K; ++i) {
            uint64_t cur = acc.at(i) + carry;
            res.at(i) = static_cast<uint32_t>(cur % q);
            carry = cur / q;
        }
        return res;
    }

    static RadixQ rq_from_T(T val, uint32_t q, size_t K)
    {
        RadixQ res(K, 0);
        for (size_t i = 0; i < K && val > T(0); ++i) {
            uint32_t rem = ::bignum::siqs::rem_u32(val, q);
            res.at(i) = rem;
            val = ::bignum::siqs::div_by_u32(val, q);
        }
        return res;
    }

    static RadixQ rq_from_signed(int64_t val, uint32_t q, size_t K)
    {
        bool neg = (val < 0);
        uint64_t abs_v = static_cast<uint64_t>(neg ? -val : val);
        RadixQ pos(K, 0);
        for (size_t i = 0; i < K && abs_v > 0; ++i) {
            pos.at(i) = static_cast<uint32_t>(abs_v % q);
            abs_v /= q;
        }
        if (!neg)
            return pos;
        RadixQ zero(K, 0);
        return rq_sub(zero, pos, q, K);
    }

    // Знаковое расширение уравновешенного представления при переходе на бóльшую
    // точность k_out > k_size: цифры старше текущей длины дополняются (q-1)
    // для отрицательного значения и 0 для положительного. Простое resize(..., 0)
    // ломало отрицательные коэффициенты: представление q^k_size - |v| перестает
    // быть congruent -v mod q^k_out.
    static RadixQ rq_extend_sign(RadixQ v, size_t k_out, uint32_t q)
    {
        if (v.size() < k_out) {
            uint32_t fill = (!v.empty() && v.back() > q / 2) ? q - 1 : 0;
            v.resize(k_out, fill);
        }
        return v;
    }

    // Умножение полиномов по модулю (f(x), q^K)
    std::vector<RadixQ> poly_mul_radix_q(const std::vector<RadixQ> &A,
                                         const std::vector<RadixQ> &B,
                                         const std::vector<RadixQ> &f_rq,
                                         uint32_t q,
                                         size_t K) const
    {
        size_t d = m_d;
        std::vector<RadixQ> prod(2 * d - 1, RadixQ(K, 0));
        for (size_t i = 0; i < d; ++i) {
            for (size_t j = 0; j < d; ++j) {
                RadixQ term = rq_mul(A.at(i), B.at(j), q, K);
                prod.at(i + j) = rq_add(prod.at(i + j), term, q, K);
            }
        }
        for (int i = static_cast<int>(2 * d - 2); i >= static_cast<int>(d); --i) {
            const RadixQ &lead = prod.at(i);
            for (size_t j = 0; j < d; ++j) {
                RadixQ sub = rq_mul(lead, f_rq.at(j), q, K);
                prod.at(i - d + j) = rq_sub(prod.at(i - d + j), sub, q, K);
            }
        }
        prod.resize(d);
        return prod;
    }

    // Быстрое дерево умножения линейных двучленов (a - b*x) по модулю (f(x), q^K).
    // Каждый уровень хранит только нужную ему точность: коэффициенты произведения
    // n линейных множителей растут как n*(бит_множителя + log2 n). Если хранить
    // все узлы сразу в полной точности K, нижние уровни дерева съедают
    // O(число_узлов * K) памяти (десятки ГБ на реальных объёмах строк).
    std::vector<RadixQ> build_product_tree_rq(const std::vector<std::pair<int64_t, uint32_t>> &pairs,
                                              const std::vector<RadixQ> &f_rq,
                                              uint32_t q,
                                              size_t K) const
    {
        size_t d = m_d;
        const size_t q_bits = static_cast<size_t>(32 - std::countl_zero(q));
        const uint64_t factor_bits = 20; // |a|, b < 2^17 с запасом

        auto level_digits = [&](size_t covered) -> size_t {
            uint64_t bits = static_cast<uint64_t>(covered)
                            * (factor_bits + (64 - std::countl_zero(covered)));
            size_t digits = (bits + q_bits - 1) / q_bits + 4;
            return std::min(K, digits);
        };

        // Знаковое расширение уравновешенного представления: см. rq_extend_sign.
        auto extend_node = [&](std::vector<RadixQ> node, size_t k_out) {
            for (auto &c : node)
                c = rq_extend_sign(std::move(c), k_out, q);
            return node;
        };

        std::vector<std::vector<RadixQ>> level;
        level.reserve(pairs.size());

        const size_t k_leaf = level_digits(1);
        for (const auto &ab : pairs) {
            std::vector<RadixQ> lin(d, RadixQ(k_leaf, 0));
            lin.at(0) = rq_from_signed(ab.first, q, k_leaf);
            lin.at(1) = rq_from_signed(-static_cast<int64_t>(ab.second), q, k_leaf);
            level.push_back(std::move(lin));
        }

        size_t covered = 1;
        while (level.size() > 1) {
            covered *= 2;
            const size_t k_lvl = level_digits(covered);
            std::vector<std::vector<RadixQ>> next_level;
            next_level.reserve((level.size() + 1) / 2);
            for (size_t i = 0; i + 1 < level.size(); i += 2) {
                auto lhs = extend_node(std::move(level.at(i)), k_lvl);
                auto rhs = extend_node(std::move(level.at(i + 1)), k_lvl);
                next_level.push_back(poly_mul_radix_q(lhs, rhs, f_rq, q, k_lvl));
            }
            if (level.size() & 1) {
                next_level.push_back(extend_node(std::move(level.back()), k_lvl));
            }
            level = std::move(next_level);
        }
        return level.at(0);
    }

    // Вычисление образа алгебраического квадратного корня через F_{q^d} и подъем Гензеля в базе q
    std::optional<T> compute_algebraic_sqrt_mod_n(const std::vector<size_t> &dep_rels) const
    {
        uint32_t q = 0;
        std::vector<uint32_t> f_mod_q(m_d + 1, 0);

        // Подбираем КРУПНЕЙШее допустимое q (< 2^16, ≡ 3 mod 4): каждая цифра
        // базы q несёт ~16 бит вместо ~11 у крошечного q≈2000, число цифр K
        // падает в ~1.5 раза, а все операции в базе q (дерево произведений,
        // подъём Ньютона) стоят O(K^2) — суммарный выигрыш ~2.5x на извлечении
        // корня. Проверка неприводимости — быстрая (Фробениус), O(d^3 log q).
        for (uint32_t cand = 65519; cand >= 2003; cand -= 4) {
            bool is_p = true;
            for (uint32_t div = 3; 1ULL * div * div <= cand; div += 2) {
                if (cand % div == 0) {
                    is_p = false;
                    break;
                }
            }
            if (!is_p)
                continue;

            for (size_t i = 0; i <= m_d; ++i) {
                f_mod_q.at(i) = ::bignum::siqs::rem_u32(m_coeffs.at(i), cand);
            }
            if (is_irreducible_mod_q_fast(f_mod_q, cand)) {
                q = cand;
                break;
            }
        }

        if (q == 0)
            return std::nullopt;

        std::vector<std::pair<int64_t, uint32_t>> pairs;
        pairs.reserve(dep_rels.size() * 2);
        for (size_t idx : dep_rels) {
            const auto &rel = m_relations.at(idx);
            pairs.emplace_back(rel.a, rel.b);
            if (rel.has_pair) {
                pairs.emplace_back(rel.a2, rel.b2);
            }
        }

        // Вычисляем производную f'(x) mod q и её квадрат f'(x)^2 mod (f(x), q)
        std::vector<uint32_t> f_deriv_q(m_d, 0);
        for (size_t i = 1; i <= m_d; ++i) {
            f_deriv_q.at(i - 1) = static_cast<uint32_t>((1ULL * i * f_mod_q.at(i)) % q);
        }
        auto f_deriv_sq_q = poly_mul_mod_q(f_deriv_q, f_deriv_q, f_mod_q, q);

        // Произведение Pi(x) = f'(x)^2 * prod(a - b*x) по модулю q
        std::vector<uint32_t> prod_q = f_deriv_sq_q;
        for (const auto &ab : pairs) {
            std::vector<uint32_t> lin(m_d, 0);
            int64_t a_m = ab.first % static_cast<int64_t>(q);
            if (a_m < 0)
                a_m += q;
            uint32_t b_m = ab.second % q;
            lin.at(0) = static_cast<uint32_t>(a_m);
            lin.at(1) = (b_m == 0) ? 0 : (q - b_m);
            prod_q = poly_mul_mod_q(prod_q, lin, f_mod_q, q);
        }

        U128 q_pow(1);
        for (uint32_t i = 0; i < m_d; ++i) {
            q_pow = q_pow * U128(q);
        }

        // --- Арифметика в поле F_{q^d} = F_q[x]/(f) ---
        auto poly_pow_field = [&](std::vector<uint32_t> base, U128 e) {
            std::vector<uint32_t> res(m_d, 0);
            res.at(0) = 1;
            while (e > U128(0)) {
                if ((e.low() & 1) != 0)
                    res = poly_mul_mod_q(res, base, f_mod_q, q);
                base = poly_mul_mod_q(base, base, f_mod_q, q);
                e >>= 1;
            }
            return res;
        };
        auto is_field_one = [&](const std::vector<uint32_t> &v) {
            if (v.at(0) != 1)
                return false;
            for (size_t i = 1; i < m_d; ++i)
                if (v.at(i) != 0)
                    return false;
            return true;
        };
        auto is_field_zero = [&](const std::vector<uint32_t> &v) {
            for (size_t i = 0; i < m_d; ++i)
                if (v.at(i) != 0)
                    return false;
            return true;
        };

        // Квадратный корень в F_{q^d}: s^2 == a. Для Q == 3 (mod 4) — простая
        // формула, иначе (всегда при чётной степени d: Q == 1 (mod 4)) —
        // алгоритм Тоннели—Шенкса. Ранее использовалась неверная степень
        // (3(Q-1)/4), дававшая s^2 == 1 вместо s^2 == a.
        auto field_sqrt = [&](const std::vector<uint32_t> &a)
            -> std::optional<std::vector<uint32_t>> {
            if (is_field_zero(a))
                return std::vector<uint32_t>(m_d, 0);

            // Q - 1 = 2^S * t, t нечётно
            U128 t = q_pow - U128(1);
            uint32_t S = 0;
            while ((t.low() & 1) == 0) {
                t >>= 1;
                ++S;
            }

            if (S == 1) { // Q == 3 (mod 4): s = a^((Q+1)/4)
                auto s = poly_pow_field(a, (q_pow + U128(1)) >> 2);
                auto sq = poly_mul_mod_q(s, s, f_mod_q, q);
                return (sq == a) ? std::optional<std::vector<uint32_t>>(s) : std::nullopt;
            }

            // Ищем квадратичный невычет z: z^((Q-1)/2) != 1
            U128 exp_legendre = t << (S - 1); // (Q-1)/2 = t * 2^(S-1)
            std::vector<uint32_t> z(m_d, 0);
            bool found = false;
            for (uint32_t c = 0; c < q && !found; ++c) {
                std::vector<uint32_t> cand(m_d, 0);
                cand.at(0) = c;
                if (m_d > 1)
                    cand.at(1) = 1; // x + c
                auto leg = poly_pow_field(cand, exp_legendre); // z^((Q-1)/2)
                if (!is_field_one(leg)) {
                    z = std::move(cand);
                    found = true;
                }
            }
            if (!found)
                return std::nullopt;

            // Тоннели—Шенкс
            std::vector<uint32_t> R = poly_pow_field(a, (t + U128(1)) >> 1); // a^((t+1)/2)
            std::vector<uint32_t> b = poly_pow_field(a, t);                  // a^t
            std::vector<uint32_t> c = poly_pow_field(z, t);                  // z^t
            uint32_t r = S;

            while (!is_field_one(b)) {
                uint32_t i = 0;
                std::vector<uint32_t> b_pow = b;
                while (i < r) {
                    if (is_field_one(b_pow))
                        break;
                    b_pow = poly_mul_mod_q(b_pow, b_pow, f_mod_q, q);
                    ++i;
                }
                if (i == 0 || i >= r)
                    return std::nullopt; // a — не квадрат

                std::vector<uint32_t> t2 = c;
                for (uint32_t k = 0; k + 1 < r - i; ++k)
                    t2 = poly_mul_mod_q(t2, t2, f_mod_q, q);
                R = poly_mul_mod_q(R, t2, f_mod_q, q);
                auto t2_sq = poly_mul_mod_q(t2, t2, f_mod_q, q);
                b = poly_mul_mod_q(b, t2_sq, f_mod_q, q);
                c = t2_sq;
                r = i;
            }

            auto sq = poly_mul_mod_q(R, R, f_mod_q, q);
            return (sq == a) ? std::optional<std::vector<uint32_t>>(R) : std::nullopt;
        };

        // inv_sqrt: s^2 == Pi^{-1}, где Pi^{-1} = Pi^{Q-2}
        auto compute_inv_sqrt = [&](const std::vector<uint32_t> &p_in)
            -> std::optional<std::vector<uint32_t>> {
            if (is_field_zero(p_in))
                return std::nullopt;
            return field_sqrt(poly_pow_field(p_in, q_pow - U128(2)));
        };

        bool negate_prod = false;
        std::vector<uint32_t> inv_sqrt_q;
        bool sqrt_field_ok = false;

        // Знак единицы в Z[alpha] заранее неизвестен: если Pi не является
        // квадратом (или дает -1), повторяем с -Pi.
        for (int attempt = 0; attempt < 2 && !sqrt_field_ok; ++attempt) {
            if (attempt == 1) {
                negate_prod = true;
                for (size_t i = 0; i < m_d; ++i) {
                    if (prod_q.at(i) != 0)
                        prod_q.at(i) = q - prod_q.at(i);
                }
            }
            auto s_opt = compute_inv_sqrt(prod_q);
            if (!s_opt)
                continue;
            auto sq_q = poly_mul_mod_q(*s_opt, *s_opt, f_mod_q, q);
            auto check_one = poly_mul_mod_q(prod_q, sq_q, f_mod_q, q);
            bool is_one = true;
            for (size_t i = 0; i < m_d; ++i) {
                if (check_one.at(i) != ((i == 0) ? 1u : 0u)) {
                    is_one = false;
                    break;
                }
            }
            if (is_one) {
                inv_sqrt_q = std::move(*s_opt);
                sqrt_field_ok = true;
            }
        }

        if (!sqrt_field_ok)
            return std::nullopt;

        // 3. Оцениваем глубину подъема K (в разрядах по основанию q).
        // Коэффициенты целочисленного корня g (если он существует) имеют порядок
        // ~2^(пара·27) бит (максимальное вложение |a - b*rho| ~ b*m); чтобы
        // центрированное представление не «заворачивалось», берём запас ×2.
        size_t total_bits_needed = (pairs.size() * (m_m.bit_width() + 20)) + m_N.bit_width() * 2
                                   + 128;
        size_t q_bits = static_cast<size_t>(std::max(10.0, std::floor(std::log2(q))));
        size_t target_K = (total_bits_needed + q_bits - 1) / q_bits + 4;

        std::vector<RadixQ> f_full(m_d + 1);
        for (size_t i = 0; i <= m_d; ++i) {
            f_full.at(i) = rq_from_T(m_coeffs.at(i), q, target_K);
        }

        // Вычисляем f'(x)^2 в базе q^target_K
        std::vector<RadixQ> f_deriv_rq(m_d, RadixQ(target_K, 0));
        for (size_t i = 1; i <= m_d; ++i) {
            RadixQ deg_rq = rq_from_signed(static_cast<int64_t>(i), q, target_K);
            f_deriv_rq.at(i - 1) = rq_mul(f_full.at(i), deg_rq, q, target_K);
        }
        auto f_deriv_sq_rq = poly_mul_radix_q(f_deriv_rq, f_deriv_rq, f_full, q, target_K);

        auto prod_full = build_product_tree_rq(pairs, f_full, q, target_K);
        for (auto &coeff : prod_full)
            coeff = rq_extend_sign(std::move(coeff), target_K, q);
        prod_full = poly_mul_radix_q(prod_full, f_deriv_sq_rq, f_full, q, target_K);

        if (negate_prod) {
            RadixQ zero_rq(target_K, 0);
            for (size_t i = 0; i < m_d; ++i) {
                prod_full.at(i) = rq_sub(zero_rq, prod_full.at(i), q, target_K);
            }
        }

        // 4. Квадратичный подъем Ньютона — Гензеля: K = 1 -> 2 -> 4 -> ... -> target_K
        size_t cur_K = 1;
        std::vector<RadixQ> S(m_d, RadixQ(1, 0));
        for (size_t i = 0; i < m_d; ++i) {
            S.at(i).at(0) = inv_sqrt_q.at(i);
        }

        uint32_t inv2_q = (q + 1) / 2;

        while (cur_K < target_K) {
            size_t next_K = std::min(cur_K * 2, target_K);

            std::vector<RadixQ> f_cur(m_d + 1, RadixQ(next_K, 0));
            for (size_t i = 0; i <= m_d; ++i) {
                std::copy_n(f_full.at(i).begin(), next_K, f_cur.at(i).begin());
            }

            std::vector<RadixQ> prod_cur(m_d, RadixQ(next_K, 0));
            for (size_t i = 0; i < m_d; ++i) {
                std::copy_n(prod_full.at(i).begin(), next_K, prod_cur.at(i).begin());
                S.at(i) = rq_extend_sign(std::move(S.at(i)), next_K, q);
            }
            auto S2 = poly_mul_radix_q(S, S, f_cur, q, next_K);
            auto P_S2 = poly_mul_radix_q(prod_cur, S2, f_cur, q, next_K);
            std::vector<RadixQ> three_minus(m_d, RadixQ(next_K, 0));
            RadixQ three_rq(next_K, 0);
            three_rq.at(0) = 3;
            RadixQ zero_rq(next_K, 0);
            RadixQ half_rq(next_K, 0);
            half_rq.at(0) = inv2_q;
            for (size_t k = 1; k < next_K; ++k) {
                half_rq.at(k) = (q - 1) / 2;
            }
            for (size_t i = 0; i < m_d; ++i) {
                if (i == 0) {
                    three_minus.at(0) = rq_sub(three_rq, P_S2.at(0), q, next_K);
                } else {
                    three_minus.at(i) = rq_sub(zero_rq, P_S2.at(i), q, next_K);
                }
                three_minus.at(i) = rq_mul(three_minus.at(i), half_rq, q, next_K);
            }
            S = poly_mul_radix_q(S, three_minus, f_cur, q, next_K);
            cur_K = next_K;
        }

        auto sqrt_rq = poly_mul_radix_q(prod_full, S, f_full, q, target_K);

#ifndef NDEBUG
        // ДИАГНОСТИКА: самопроверка подъёма — sqrt_rq^2 == P (mod q^K) поцифрово.
        // Если сломан подъём — BROKEN; если OK, а Y^2 != X^2, то единица в Z[alpha]
        // не является квадратом (подъём даёт q-адический, а не целочисленный корень)
        // либо целочисленный корень не поместился в q^K (тогда used_bits ~= K*q_bits).
        {
            auto sq_full = poly_mul_radix_q(sqrt_rq, sqrt_rq, f_full, q, target_K);
            bool ok_K = true;
            for (size_t i = 0; i < m_d && ok_K; ++i) {
                if (sq_full.at(i) != prod_full.at(i))
                    ok_K = false;
            }
            // «Использованные биты» центрированного корня: сколько старших цифр
            // занято знаковыми разрядами (0 — для положительных, q-1 — для отрицательных).
            size_t q_bits_diag = static_cast<size_t>(std::max(10.0, std::floor(std::log2(q))));
            size_t used_bits = 0;
            for (size_t i = 0; i < m_d; ++i) {
                const RadixQ &c = sqrt_rq.at(i);
                uint32_t fill = (c.back() > q / 2) ? (q - 1) : 0;
                size_t first = 0;
                while (first < target_K && c.at(target_K - 1 - first) == fill)
                    ++first;
                size_t bits_here = (target_K - first) * q_bits_diag;
                if (bits_here > used_bits)
                    used_bits = bits_here;
            }
            std::fprintf(stderr,
                         "DIAG sqrt: q=%u target_K=%zu (%zu бит) pairs=%zu negate=%d "
                         "sqrt^2==P(mod q^K):%s used_bits=%zu\n",
                         q, target_K, target_K * q_bits_diag, pairs.size(),
                         static_cast<int>(negate_prod), ok_K ? "OK" : "BROKEN", used_bits);
            std::fflush(stderr);
        }
#endif // NDEBUG

        // 5. Переводим центрированные коэффициенты из базы q в Z_N и подставляем x = m
        T Y = T(0);
        T m_pow = T(1);
        for (size_t i = 0; i < m_d; ++i) {
            const RadixQ &c_rq = sqrt_rq.at(i);
            bool is_neg_coeff = (c_rq.back() > q / 2);
            RadixQ abs_rq = is_neg_coeff ? rq_sub(RadixQ(target_K, 0), c_rq, q, target_K) : c_rq;

            T c_mod_n = T(0);
            for (int k = static_cast<int>(target_K) - 1; k >= 0; --k) {
                ::bignum::siqs::mult_mod(c_mod_n, T(q), m_N);
                c_mod_n = c_mod_n + T(abs_rq.at(k));
                if (c_mod_n >= m_N)
                    c_mod_n = c_mod_n - m_N;
            }
            ::bignum::siqs::mult_mod(c_mod_n, m_pow, m_N);
            if (!is_neg_coeff) {
                Y = Y + c_mod_n;
                if (Y >= m_N)
                    Y = Y - m_N;
            } else {
                Y = (Y >= c_mod_n) ? (Y - c_mod_n) : (m_N - (c_mod_n - Y));
            }
            ::bignum::siqs::mult_mod(m_pow, m_m, m_N);
        }
        return Y;
    }

    std::optional<T> solve_linear_algebra(std::function<void(const std::string &)> status_cb)
    {
        size_t rfb_len = m_RFB.size();
        size_t afb_len = m_AFB.size();
        size_t qcb_len = m_QCB.size();

        size_t base_prime_cols = 1 + rfb_len + 1 + afb_len;
        size_t total_cols = base_prime_cols + qcb_len;
        size_t M = m_relations.size();

        if (M <= total_cols) {
            if (status_cb)
                status_cb("GNFS: Недостаточно соотношений для построения матрицы ("
                          + std::to_string(M) + " <= " + std::to_string(total_cols) + ")");
            return std::nullopt;
        }

        if (status_cb)
            status_cb("GNFS: Построение бинарной матрицы: " + std::to_string(M) + " строк x "
                      + std::to_string(total_cols) + " столбцов");

        size_t words = (total_cols + 63) / 64;
        size_t hist_words = (M + 63) / 64;

        struct MatrixRow
        {
            std::vector<uint64_t> mask;
            std::vector<uint64_t> history;
        };

        std::vector<MatrixRow> matrix(M);
        size_t alg_offset = 1 + rfb_len;
        size_t qcb_offset = base_prime_cols;

        for (size_t r = 0; r < M; ++r) {
            const auto &rel = m_relations[r];
            matrix[r].mask.assign(words, 0);
            matrix[r].history.assign(hist_words, 0);
            matrix[r].history[r / 64] |= (1ULL << (r % 64));

            auto toggle_col = [&](size_t col) { matrix[r].mask[col / 64] ^= (1ULL << (col % 64)); };

            if (rel.rat_neg)
                toggle_col(0);
            for (size_t i = 0; i < rel.rat_idx.size(); ++i) {
                if (rel.rat_pow[i] & 1)
                    toggle_col(1 + rel.rat_idx[i]);
            }

            if (rel.alg_neg)
                toggle_col(alg_offset);
            for (size_t i = 0; i < rel.alg_idx.size(); ++i) {
                if (rel.alg_pow[i] & 1)
                    toggle_col(alg_offset + 1 + rel.alg_idx[i]);
            }

            for (size_t i = 0; i < qcb_len; ++i) {
                uint32_t q = m_QCB[i].q;
                uint32_t s = m_QCB[i].s;

                auto eval_qcb = [&](int64_t a_val, uint32_t b_val) {
                    int64_t a_mod = a_val % static_cast<int64_t>(q);
                    if (a_mod < 0)
                        a_mod += q;
                    uint64_t bs_mod = (1ULL * b_val * s) % q;
                    uint32_t v = static_cast<uint32_t>((static_cast<uint64_t>(a_mod) + q - bs_mod)
                                                       % q);
                    return ::bignum::siqs::legendre_symbol(v, q) == -1;
                };

                bool qcb_bit = eval_qcb(rel.a, rel.b);
                if (rel.has_pair) {
                    qcb_bit ^= eval_qcb(rel.a2, rel.b2);
                }
                if (qcb_bit)
                    toggle_col(qcb_offset + i);
            }
        }

        // Прямой и обратный ход метода Гаусса над GF(2)
        size_t row = 0;
        for (size_t col = 0; col < total_cols && row < M; ++col) {
            size_t col_w = col / 64;
            uint64_t col_mask = 1ULL << (col % 64);
            size_t pivot = row;
            while (pivot < M && !(matrix[pivot].mask[col_w] & col_mask)) {
                pivot++;
            }
            if (pivot >= M)
                continue;
            if (row != pivot)
                std::swap(matrix[row], matrix[pivot]);

            const auto &pivot_mask = matrix[row].mask;
            const auto &pivot_hist = matrix[row].history;
            for (size_t r = 0; r < M; ++r) {
                if (r != row && (matrix[r].mask[col_w] & col_mask)) {
                    for (size_t w = col_w; w < words; ++w)
                        matrix[r].mask[w] ^= pivot_mask[w];
                    for (size_t hw = 0; hw < hist_words; ++hw)
                        matrix[r].history[hw] ^= pivot_hist[hw];
                }
            }
            row++;
        }

        if (status_cb) {
            status_cb("GNFS: СЛАУ: Гаусс завершён, перебор векторов ядра...");
        }

        size_t null_rows = 0;
        size_t sqrt_ok = 0;
        size_t sq_match = 0;

        // Нулевые строки = базис ядра. Базисные векторы МОГУТ все давать один и
        // тот же «неквадратный» класс (единица/класс кольца Z[alpha]) — тогда ни
        // один прямой вектор не даст целочисленного корня, хотя комбинации дают.
        // Любая сумма векторов ядра тоже вектор ядра, поэтому перебираем и их:
        // v0, v1, v0^v1, затем для j>=2: vj, vj^v0, vj^v1, vj^v0^v1 — это покрывает
        // все 4 сочетания по первым двум направлениям ядра для каждого j.
        std::vector<size_t> null_list;
        for (size_t r = 0; r < M; ++r) {
            bool all_zero = true;
            for (size_t w = 0; w < words; ++w) {
                if (matrix[r].mask[w] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero)
                null_list.push_back(r);
        }

        std::vector<uint64_t> hist_comb(hist_words, 0);
        const size_t total_attempts =
            3 + (null_list.size() > 2 ? (null_list.size() - 2) * 4 : 0);
        for (size_t attempt = 0; attempt < total_attempts; ++attempt) {
            if (::u128::Globals::LoadStop())
                break;
            size_t j, mask;
            if (attempt == 0) {
                j = 0;
                mask = 0; // v0
            } else if (attempt == 1) {
                j = 1;
                mask = 0; // v1
            } else if (attempt == 2) {
                j = 0;
                mask = 2; // v0 ^ v1
            } else {
                j = (attempt - 3) / 4 + 2;
                mask = (attempt - 3) % 4;
            }
            if (j >= null_list.size())
                break;
            if ((mask & 3) != 0 && null_list.size() < 2)
                break;

            null_rows++;
            // comb = v_j ^ (mask&1 ? v0 : 0) ^ (mask&2 ? v1 : 0)
            std::copy(matrix[null_list[j]].history.begin(),
                      matrix[null_list[j]].history.end(), hist_comb.begin());
            if (mask & 1) {
                for (size_t w = 0; w < hist_words; ++w)
                    hist_comb[w] ^= matrix[null_list[0]].history[w];
            }
            if (mask & 2) {
                for (size_t w = 0; w < hist_words; ++w)
                    hist_comb[w] ^= matrix[null_list[1]].history[w];
            }

            std::vector<uint64_t> rfb_powers(rfb_len, 0);
            std::unordered_map<uint64_t, uint64_t> lp_rat_counts;
            std::vector<size_t> dep_rels;

            for (size_t local_idx = 0; local_idx < M; ++local_idx) {
                if (hist_comb[local_idx / 64] & (1ULL << (local_idx % 64))) {
                    dep_rels.push_back(local_idx);
                    const auto &rel = m_relations[local_idx];
                    for (size_t i = 0; i < rel.rat_idx.size(); ++i) {
                        rfb_powers[rel.rat_idx[i]] += rel.rat_pow[i];
                    }
                    if (rel.has_pair && rel.large_rat_p > 1) {
                        lp_rat_counts[rel.large_rat_p] += 2;
                    }
                }
            }
            if (dep_rels.empty())
                continue;

            T X = T(1);
            for (size_t i = 0; i < rfb_len; ++i) {
                uint64_t half_p = rfb_powers[i] / 2;
                if (half_p > 0) {
                    T p_val(m_RFB[i].p);
                    ::bignum::siqs::int_power_mod(p_val, T(half_p), m_N);
                    ::bignum::siqs::mult_mod(X, p_val, m_N);
                }
            }
            for (const auto &[lp, cnt] : lp_rat_counts) {
                uint64_t half_lp = cnt / 2;
                if (half_lp > 0) {
                    T p_val(lp);
                    ::bignum::siqs::int_power_mod(p_val, T(half_lp), m_N);
                    ::bignum::siqs::mult_mod(X, p_val, m_N);
                }
            }

            T f_deriv_m = T(0);
            T m_pow = T(1);
            for (size_t i = 1; i <= m_d; ++i) {
                T term = ::bignum::siqs::div_r(m_coeffs.at(i), m_N);
                ::bignum::siqs::mult_mod(term, T(static_cast<uint64_t>(i)), m_N);
                ::bignum::siqs::mult_mod(term, m_pow, m_N);
                f_deriv_m = f_deriv_m + term;
                if (f_deriv_m >= m_N)
                    f_deriv_m = f_deriv_m - m_N;
                ::bignum::siqs::mult_mod(m_pow, m_m, m_N);
            }
            ::bignum::siqs::mult_mod(X, f_deriv_m, m_N);

            if (status_cb && (null_rows <= 3 || null_rows % 25 == 0)) {
                status_cb("GNFS: СЛАУ: вектор ядра #" + std::to_string(null_rows)
                          + " (строк в связи: " + std::to_string(dep_rels.size())
                          + "), извлечение корня...");
            }

            auto y_opt = compute_algebraic_sqrt_mod_n(dep_rels);
            if (!y_opt.has_value()) {
                if (status_cb && (null_rows <= 3 || null_rows % 25 == 0))
                    status_cb("GNFS: СЛАУ: вектор ядра #" + std::to_string(null_rows)
                              + ": корень не извлечён (sqrt=" + std::to_string(sqrt_ok) + ")");
                continue;
            }
            sqrt_ok++;
            T Y = *y_opt;

            T X2 = X;
            ::bignum::siqs::square_mod(X2, m_N);
            T Y2 = Y;
            ::bignum::siqs::square_mod(Y2, m_N);

#ifndef NDEBUG
            // ДИАГНОСТИКА: прямой пересчёт P(m) = f'(m)^2 * prod(a - b*m) (mod N).
            // Сверка X^2 и Y^2 с ним локализует ошибку: X^2!=P(m) — сторона X/матрицы,
            // Y^2!=P(m) — извлечение корня.
            if (status_cb && null_rows <= 3) {
                // Точная сверка факторизаций соотношений в связи (без переполнений —
                // делением): ∏p^e · LP^2 == |F(a,b)| (·|F(a2,b2)|). Несовпадение
                // означает баг разложения/спаривания (элементарная ошибка), совпадение
                // при Y^2!=P(m) — «неквадратная» единица в Z[alpha] или баг подъёма.
                size_t bad_rat = 0, bad_alg = 0;
                auto divide_out = [](T &q, uint32_t p, uint16_t &e) {
                    while (e > 0 && ::bignum::siqs::rem_u32(q, p) == 0) {
                        q = ::bignum::siqs::div_by_u32(q, p);
                        --e;
                    }
                };
                for (size_t idx : dep_rels) {
                    const auto &rel2 = m_relations[idx];
                    // Алгебраическая сторона: |F(a,b)| и |F(a2,b2)| против alg_idx/alg_pow + LP^2
                    bool f_neg = false;
                    T fa = eval_homogeneous_poly(m_coeffs, rel2.a, rel2.b, f_neg);
                    T fb(1);
                    if (rel2.has_pair)
                        fb = eval_homogeneous_poly(m_coeffs, rel2.a2, rel2.b2, f_neg);
                    bool ok_a = true;
                    for (size_t k = 0; k < rel2.alg_idx.size() && ok_a; ++k) {
                        uint16_t e = rel2.alg_pow[k];
                        uint32_t p = m_AFB.at(rel2.alg_idx[k]).p;
                        divide_out(fa, p, e);
                        if (e > 0 && fb != T(1)) {
                            divide_out(fb, p, e);
                        }
                        if (e > 0)
                            ok_a = false;
                    }
                    if (ok_a && rel2.large_alg_p > 1) {
                        uint32_t p32 = static_cast<uint32_t>(rel2.large_alg_p);
                        uint16_t e2 = 2;
                        divide_out(fa, p32, e2);
                        if (e2 > 0 && fb != T(1))
                            divide_out(fb, p32, e2);
                        if (e2 > 0)
                            ok_a = false;
                    }
                    if (ok_a && (fa != T(1) || fb != T(1)))
                        ok_a = false;
                    if (!ok_a)
                        ++bad_alg;

                    // Рациональная сторона: ∏|a-bm| против rat_idx/rat_pow + LP^2
                    T ra(static_cast<uint64_t>(rel2.a >= 0 ? rel2.a : -rel2.a));
                    T rbv(T(rel2.b) * m_m);
                    T rat_abs = (ra >= rbv) ? (ra - rbv) : (rbv - ra);
                    T rat2(1);
                    if (rel2.has_pair) {
                        T ra2(static_cast<uint64_t>(rel2.a2 >= 0 ? rel2.a2 : -rel2.a2));
                        T rb2v(T(rel2.b2) * m_m);
                        rat2 = (ra2 >= rb2v) ? (ra2 - rb2v) : (rb2v - ra2);
                    }
                    bool ok_r = true;
                    for (size_t k = 0; k < rel2.rat_idx.size() && ok_r; ++k) {
                        uint16_t e = rel2.rat_pow[k];
                        uint32_t p = m_RFB.at(rel2.rat_idx[k]).p;
                        divide_out(rat_abs, p, e);
                        if (e > 0 && rat2 != T(1))
                            divide_out(rat2, p, e);
                        if (e > 0)
                            ok_r = false;
                    }
                    if (ok_r && rel2.large_rat_p > 1) {
                        uint32_t p32 = static_cast<uint32_t>(rel2.large_rat_p);
                        uint16_t e2 = 2;
                        divide_out(rat_abs, p32, e2);
                        if (e2 > 0 && rat2 != T(1))
                            divide_out(rat2, p32, e2);
                        if (e2 > 0)
                            ok_r = false;
                    }
                    if (ok_r && (rat_abs != T(1) || rat2 != T(1)))
                        ok_r = false;
                    if (!ok_r)
                        ++bad_rat;
                }

                size_t pair_cnt = 0;
                size_t neg_cnt = 0;
                T prod_m(1);
                for (size_t idx : dep_rels) {
                    const auto &rel2 = m_relations[idx];
                    auto mul_ab = [&](int64_t aa, uint32_t bb) {
                        ++pair_cnt;
                        T bmv = T(bb) * m_m;
                        T val;
                        if (aa >= 0) {
                            T av(static_cast<uint64_t>(aa));
                            if (av >= bmv) {
                                val = av - bmv;
                            } else {
                                val = m_N - (bmv - av);
                                ++neg_cnt;
                            }
                        } else {
                            val = m_N - (bmv + T(static_cast<uint64_t>(-aa)));
                            ++neg_cnt;
                        }
                        ::bignum::siqs::mult_mod(prod_m, val, m_N);
                    };
                    mul_ab(rel2.a, rel2.b);
                    if (rel2.has_pair)
                        mul_ab(rel2.a2, rel2.b2);
                }
                T fd2 = f_deriv_m;
                ::bignum::siqs::square_mod(fd2, m_N);
                ::bignum::siqs::mult_mod(prod_m, fd2, m_N);
                bool x_ok = (X2 == prod_m);
                bool y_ok = (Y2 == prod_m);
                T neg_prod_m = (prod_m == T(0)) ? T(0) : (m_N - prod_m); // -P(m) mod N
                bool y_neg_ok = (Y2 == neg_prod_m);
                status_cb("GNFS: ДИАГ: вектор #" + std::to_string(null_rows)
                          + " пар=" + std::to_string(pair_cnt)
                          + " отр.норм=" + std::to_string(neg_cnt)
                          + " битых разл. рац/алг=" + std::to_string(bad_rat) + "/" + std::to_string(bad_alg)
                          + " X^2==P(m):" + (x_ok ? "ДА" : "НЕТ")
                          + " Y^2==P(m):" + (y_ok ? "ДА" : "НЕТ")
                          + " Y^2==-P(m):" + (y_neg_ok ? "ДА" : "НЕТ"));
            }
#endif // NDEBUG

            if (X2 == Y2)
                sq_match++;

            if (status_cb && (null_rows <= 5 || null_rows % 25 == 0)) {
                status_cb("GNFS: СЛАУ: вектор ядра #" + std::to_string(null_rows)
                          + ": X^2==Y^2 " + (X2 == Y2 ? "ДА" : "НЕТ") + " (match="
                          + std::to_string(sq_match) + "/" + std::to_string(sqrt_ok) + ")");
            }

            T diff = (X >= Y) ? (X - Y) : (Y - X);
            T sum = ::bignum::siqs::div_r(X + Y, m_N);
            T g1 = ::bignum::siqs::gcd_generic(diff, m_N);
            if (g1 > T(1) && g1 < m_N) {
                if (status_cb)
                    status_cb("GNFS: Успех! Найден делитель: " + g1.toString());
                return g1;
            }
            T g2 = ::bignum::siqs::gcd_generic(sum, m_N);
            if (g2 > T(1) && g2 < m_N) {
                if (status_cb)
                    status_cb("GNFS: Успех! Найден делитель: " + g2.toString());
                return g2;
            }
        }

#ifndef NDEBUG
        if (status_cb) {
            status_cb("DEBUG: rels=" + std::to_string(M) + " | cols=" + std::to_string(total_cols)
                      + " | null=" + std::to_string(null_rows) + " | sqrt="
                      + std::to_string(sqrt_ok) + " | match=" + std::to_string(sq_match));
            std::this_thread::sleep_for(std::chrono::seconds(6));
        }
#endif // NDEBUG
        return std::nullopt;
    }
};
template<typename T>
inline std::optional<T> number_field_sieve(
    const T &n, std::function<void(const std::string &)> status_cb = nullptr)
{
    GNFSSolver<T> solver(n);
    return solver.solve(status_cb);
}
} // namespace bignum::gnfs
