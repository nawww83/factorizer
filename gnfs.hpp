#pragma once

#include "siqs.hpp"
#include "u128.hpp"
#include "u128_utils.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
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

//
// --- Основной класс Решета Числового Поля (GNFS) ---

template<typename T>
class GNFSSolver
{
    friend struct GNFSProbe; // Доступ для диагностического пробника/тестов

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

    std::mutex m_rel_mutex;
    std::vector<GNFSRelation<T>> m_relations;
    std::unordered_set<uint64_t> m_seen_ab;

    // Таблицы для спаривания одиночных больших простых:
    std::unordered_map<uint64_t, GNFSRelation<T>> m_rat_partials;
    std::unordered_map<uint64_t, GNFSRelation<T>> m_alg_partials;

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

        // Для 128-180 бит степень d=4 дает минимальные нормы и взрывную плотность:
        if (bits < 180)
            m_d = 4;
        else if (bits < 280)
            m_d = 5;
        else
            m_d = 5;

        if (status_cb)
            status_cb("GNFS: Подбор Base-m полинома (d=" + std::to_string(m_d) + ")...");
        if (!select_polynomial()) {
            return std::nullopt;
        }

        // Подбор баз эмпирический (см. замеры для 128-бит): слишком маленькие базы
        // (800/850) дают вероятность гладкости алгебраической нормы ~0.1%,
        // и решето физически не успевает набрать строки под СЛАУ.
        // 3000/5000: ~1.4x запаса строк к b ~ 1600-2600.
        size_t rfb_target = 3000;
        size_t afb_target = 5000;
        // Полезная область малых |a| (< 65536) и малых b: алгебраическая норма
        // растёт как b^4, дальше решето даёт лишь мусор.
        m_sieve_A = 65536;

        if (status_cb)
            status_cb("GNFS: Построение рациональной и алгебраической баз...");
        if (!build_factor_bases(rfb_target, afb_target, 16)) {
            return std::nullopt;
        }

        size_t total_cols = (1 + m_RFB.size() + 1 + m_AFB.size()) + m_QCB.size();
        // Строк должно быть заметно больше числа столбцов, иначе в СЛАУ не будет
        // нетривиального вектора в ядре.
        m_relations_needed = total_cols + std::max<size_t>(64, total_cols / 64);

        m_relations.clear();
        m_seen_ab.clear();
        m_rat_partials.clear();
        m_alg_partials.clear();
        m_full_relations_count = 0;
        m_combined_count = 0;
        m_done.store(false);

        if (status_cb)
            status_cb("GNFS: Запуск двухстороннего решета (RFB + AFB + Large Primes)...");
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
    bool select_polynomial()
    {
        // Подбираем m так, чтобы m^d <= N < 2 * m^d, тогда старший коэффициент c_d РОВНО равен 1!
        m_m = nth_root(m_N, m_d);
        if (m_m <= T(2))
            return false;

        while (true) {
            T m_pow_d = T(1);
            for (uint32_t i = 0; i < m_d; ++i)
                m_pow_d = m_pow_d * m_m;

            auto [top_q, top_r] = ::bignum::siqs::div_rem(m_N, m_pow_d);
            if (top_q == T(1))
                break; // Идеально: c_d == 1 и все остальные c_i < m!
            if (top_q == T(0)) {
                m_m = m_m - T(1);
            } else {
                // Увеличиваем m, пока c_d не станет равным 1
                T step = m_m >> 4;
                if (step == T(0))
                    step = T(1);
                m_m = m_m + step;
            }
        }

        m_coeffs.clear();
        T temp = m_N;
        for (uint32_t i = 0; i < m_d; ++i) {
            auto [q, r] = ::bignum::siqs::div_rem(temp, m_m);
            m_coeffs.push_back(r);
            temp = q;
        }
        m_coeffs.push_back(temp);

        return (temp == T(1));
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

        // Не раздуваем пул LP, чтобы они постоянно сталкивались в пары:
        m_L_rat_max = max_rfb_p * 250ULL; // ~1.6 млн
        m_L_alg_max = max_afb_p * 250ULL; // ~1.8 млн

        return true;
    }

    void sieve_worker(std::atomic<uint32_t> &b_counter,
                      std::function<void(const std::string &)> status_cb)
    {
        const int M = m_sieve_A;
        std::vector<uint8_t> rat_sieve(2 * M, 0);
        std::vector<uint8_t> alg_sieve(2 * M, 0);

        // Разрешаем остаткам после решета быть больше,
        // давая возможность находить числа с одним Large Prime (до p_max^2):
        uint8_t l_rat_slack = m_RFB.back().log_p * 3 + 8;
        uint8_t l_alg_slack = m_AFB.back().log_p * 3 + 12;

        constexpr uint32_t MAX_B = 60000;
        while (!m_done.load(std::memory_order_relaxed) && !::u128::Globals::LoadStop()) {
            uint32_t b = b_counter.fetch_add(1);
            if (b > MAX_B)
                break;

            if (status_cb && (b % 10 == 0)) {
                std::lock_guard<std::mutex> lock(m_rel_mutex);
                if (!m_done.load(std::memory_order_relaxed)) {
                    status_cb("GNFS: Сбор соотношений: " + std::to_string(m_relations.size())
                              + " / " + std::to_string(m_relations_needed)
                              + " (пар: " + std::to_string(m_combined_count)
                              + ", b=" + std::to_string(b) + ")");
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

            // Пороги подобраны под базы ~3000/5000 (см. замеры): заниженные
            // пороги (22/27) пропускали в перебор на деление ~10x больше кандидатов,
            // почти не добавляя строк. Гладкие соотношения проходят с большим запасом.
            const uint8_t rat_thresh = 28;
            const uint8_t alg_thresh = 34;

            for (int blk = 0; blk < NUM_BLOCKS; ++blk) {
                int j_start = blk * block_len;
                int j_end = (blk == NUM_BLOCKS - 1) ? (2 * M) : (j_start + block_len);

                static std::atomic<uint64_t> stat_sieve_hits{0};
                static std::atomic<uint64_t> stat_rfb_fail{0};
                static std::atomic<uint64_t> stat_afb_fail{0};
                static std::atomic<uint64_t> stat_success{0};

                for (int j = j_start; j < j_end; ++j) {
                    if (rat_sieve.at(j) >= rat_thresh && alg_sieve.at(j) >= alg_thresh) {
                        stat_sieve_hits++;
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

                        std::vector<uint32_t> r_idx;
                        std::vector<uint16_t> r_pow;
                        uint64_t large_rat = 1;
                        if (!factor_over_rfb(rat_norm, r_idx, r_pow, large_rat)) {
                            stat_rfb_fail++;
                            continue;
                        }

                        bool alg_neg = false;
                        T alg_norm = eval_homogeneous_poly(m_coeffs, a, b, alg_neg);
                        if (alg_norm == T(0))
                            continue;

                        std::vector<uint32_t> a_idx;
                        std::vector<uint16_t> a_pow;
                        uint64_t large_alg = 1;
                        uint32_t large_alg_r = 0;
                        if (!factor_over_afb(alg_norm, a, b, a_idx, a_pow, large_alg, large_alg_r)) {
                            stat_afb_fail++;
                            continue;
                        }

                        // Отбрасываем соотношения с двумя большими простыми одновременно:
                        if (large_rat > 1 && large_alg > 1)
                            continue;

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

                            // 1. Полностью гладкое соотношение (0-LP)
                            if (large_rat == 1 && large_alg == 1) {
                                m_relations.push_back(std::move(rel));
                                m_full_relations_count++;
                                stat_success++;
                            }
                            // 2. Рациональный Large Prime (1-LP)
                            else if (large_rat > 1) {
                                auto it = m_rat_partials.find(large_rat);
                                if (it == m_rat_partials.end()) {
                                    m_rat_partials.emplace(large_rat, std::move(rel));
                                } else {
                                    // Найдена пара! Объединяем их в одну выживающую строку:
                                    GNFSRelation<T> combined;
                                    combined.a = it->second.a;
                                    combined.b = it->second.b;
                                    combined.has_pair = true;
                                    combined.a2 = rel.a;
                                    combined.b2 = rel.b;
                                    combined.rat_neg = (it->second.rat_neg != rel.rat_neg);
                                    combined.alg_neg = (it->second.alg_neg != rel.alg_neg);

                                    merge_factor_lists(it->second.rat_idx,
                                                       it->second.rat_pow,
                                                       rel.rat_idx,
                                                       rel.rat_pow,
                                                       combined.rat_idx,
                                                       combined.rat_pow);
                                    merge_factor_lists(it->second.alg_idx,
                                                       it->second.alg_pow,
                                                       rel.alg_idx,
                                                       rel.alg_pow,
                                                       combined.alg_idx,
                                                       combined.alg_pow);

                                    combined.large_rat_p
                                        = large_rat; // Квадрат L^2 взаимно уничтожится mod 2
                                    combined.large_alg_p = 1;

                                    m_relations.push_back(std::move(combined));
                                    m_rat_partials.erase(it);
                                    m_combined_count++;
                                    stat_success++;
                                }
                            }
                            // 3. Алгебраический Large Prime (1-LP)
                            else if (large_alg > 1) {
                                uint64_t alg_key = (large_alg << 32) | large_alg_r;
                                auto it = m_alg_partials.find(alg_key);
                                if (it == m_alg_partials.end()) {
                                    m_alg_partials.emplace(alg_key, std::move(rel));
                                } else {
                                    GNFSRelation<T> combined;
                                    combined.a = it->second.a;
                                    combined.b = it->second.b;
                                    combined.has_pair = true;
                                    combined.a2 = rel.a;
                                    combined.b2 = rel.b;
                                    combined.rat_neg = (it->second.rat_neg != rel.rat_neg);
                                    combined.alg_neg = (it->second.alg_neg != rel.alg_neg);

                                    merge_factor_lists(it->second.rat_idx,
                                                       it->second.rat_pow,
                                                       rel.rat_idx,
                                                       rel.rat_pow,
                                                       combined.rat_idx,
                                                       combined.rat_pow);
                                    merge_factor_lists(it->second.alg_idx,
                                                       it->second.alg_pow,
                                                       rel.alg_idx,
                                                       rel.alg_pow,
                                                       combined.alg_idx,
                                                       combined.alg_pow);

                                    combined.large_rat_p = 1;
                                    combined.large_alg_p = large_alg;
                                    combined.large_alg_r = large_alg_r;

                                    m_relations.push_back(std::move(combined));
                                    m_alg_partials.erase(it);
                                    m_combined_count++;
                                }
                            }

                            size_t total_ready = m_full_relations_count + m_combined_count;
                            if (status_cb && total_ready % 20 == 0) {
                                status_cb("GNFS: Сбор соотношений: " + std::to_string(total_ready)
                                          + " / " + std::to_string(m_relations_needed)
                                          + " (чистых: " + std::to_string(m_full_relations_count)
                                          + ", пар: " + std::to_string(m_combined_count) + ")");
                            }

                            // Останавливаемся ТОЛЬКО когда готовых (гарантированно парных) строк достаточно!
                            if (total_ready >= m_relations_needed) {
                                m_done.store(true, std::memory_order_relaxed);
                                return;
                            }
                        }
                    }
                }
            }
        }
        if (b_counter.load() > MAX_B && !m_done.load()) {
            if (status_cb) {
                status_cb("GNFS: Исчерпан лимит b. Собрано строк: "
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

        for (uint32_t cand = 2003; cand < 50000; cand += 4) {
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
            if (is_irreducible_mod_q(f_mod_q, cand)) {
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

        // 3. Оцениваем глубину подъема K (в разрядах по основанию q)
        size_t total_bits_needed = (pairs.size() * (m_m.bit_width() + 20)) / 2 + m_N.bit_width() * 2
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

        for (size_t r = 0; r < M; ++r) {
            bool all_zero = true;
            for (size_t w = 0; w < words; ++w) {
                if (matrix[r].mask[w] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero)
                continue;

            null_rows++;
            std::vector<uint64_t> rfb_powers(rfb_len, 0);
            std::unordered_map<uint64_t, uint64_t> lp_rat_counts;
            std::vector<size_t> dep_rels;

            for (size_t local_idx = 0; local_idx < M; ++local_idx) {
                if (matrix[r].history[local_idx / 64] & (1ULL << (local_idx % 64))) {
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

        if (status_cb) {
            status_cb("DEBUG: rels=" + std::to_string(M) + " | cols=" + std::to_string(total_cols)
                      + " | null=" + std::to_string(null_rows) + " | sqrt="
                      + std::to_string(sqrt_ok) + " | match=" + std::to_string(sq_match));
            std::this_thread::sleep_for(std::chrono::seconds(6));
        }
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
