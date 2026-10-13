#pragma once

#include "u128.hpp"
#include "ubig.hpp"
#include "u128_utils.h"
#include <random>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <optional>
#include <thread>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <iterator>
#include <string>

namespace bignum::siqs {

// --- Универсальное извлечение младших бит из любой иерархии UBig ---

template <typename U>
inline constexpr uint64_t to_u64_low(const U& val) {
    if constexpr (std::is_integral_v<U>) {
        return static_cast<uint64_t>(val);
    } else if constexpr (requires { val.low(); }) {
        return to_u64_low(val.low());
    } else {
        return 0;
    }
}

template <typename U>
inline constexpr uint32_t to_u32_low(const U& val) {
    return static_cast<uint32_t>(to_u64_low(val));
}

// --- Извлечение половинок и сборка больших чисел ---

template <typename U>
inline constexpr auto get_hi(const U& val) {
    if constexpr (requires { val.high(); }) {
        return val.high();
    } else if constexpr (requires { val.hi(); }) {
        return val.hi();
    } else if constexpr (std::is_same_v<U, ::bignum::u128::U128>) {
        return to_u64_low(val >> 64);
    }
}

template <typename U>
inline constexpr auto get_lo(const U& val) {
    if constexpr (requires { val.low(); }) {
        return val.low();
    } else if constexpr (requires { val.lo(); }) {
        return val.lo();
    } else if constexpr (std::is_same_v<U, ::bignum::u128::U128>) {
        return to_u64_low(val);
    }
}

template <typename T, typename Half>
inline T make_big(const Half& hi, const Half& lo) {
    if constexpr (std::is_same_v<T, ::bignum::u128::U128>) {
        static const bool hi_first = []() {
            ::bignum::u128::U128 test_obj{1ULL, 2ULL};
            return (to_u64_low(test_obj) == 2ULL);
        }();
        if (hi_first) {
            return ::bignum::u128::U128{hi, lo};
        } else {
            return ::bignum::u128::U128{lo, hi};
        }
    } else {
        static const bool hi_first = []() {
            Half test_h(1), test_l(2);
            T test_obj(test_h, test_l);
            return (get_lo(test_obj) == test_l);
        }();
        if (hi_first) {
            return T(hi, lo);
        } else {
            return T(lo, hi);
        }
    }
}

// Быстрый хэшер больших чисел без строковой конверсии
template <typename T>
struct BigHash {
    size_t operator()(const T& val) const noexcept {
        uint64_t lo = to_u64_low(val);
        uint64_t hi = to_u64_low(get_hi(val));
        uint64_t seed = lo ^ (hi + 0x9e3779b97f4a7c15ULL + (lo << 6) + (lo >> 2));
        seed ^= (seed >> 30) * 0xbf58476d1ce4e5b9ULL;
        seed ^= (seed >> 27) * 0x94d049bb133111ebULL;
        seed ^= (seed >> 31);
        return static_cast<size_t>(seed);
    }
};

inline uint32_t ctz32(uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<uint32_t>(__builtin_ctz(x));
#else
    uint32_t c = 0;
    while ((x & 1) == 0) { x >>= 1; c++; }
    return c;
#endif
}

inline uint32_t rem_u32_step(uint64_t val, uint32_t carry, uint32_t d) {
    if (carry == 0) {
        return static_cast<uint32_t>(val % d);
    }
#if defined(__SIZEOF_INT128__)
    unsigned __int128 cur = (static_cast<unsigned __int128>(carry) << 64) | val;
    return static_cast<uint32_t>(cur % d);
#else
    uint32_t hi = static_cast<uint32_t>(val >> 32);
    uint32_t lo = static_cast<uint32_t>(val);
    uint64_t cur1 = (static_cast<uint64_t>(carry) << 32) | hi;
    carry = static_cast<uint32_t>(cur1 % d);
    uint64_t cur2 = (static_cast<uint64_t>(carry) << 32) | lo;
    return static_cast<uint32_t>(cur2 % d);
#endif
}

template <typename T>
inline uint32_t rem_u32_step_generic(const T& val, uint32_t carry, uint32_t d) {
    if constexpr (std::is_integral_v<T>) {
        return rem_u32_step(static_cast<uint64_t>(val), carry, d);
    } else if constexpr (std::is_same_v<T, ::bignum::u128::U128>) {
        uint64_t hi = get_hi(val);
        uint64_t lo = get_lo(val);

        if (carry == 0 && hi == 0) {
            return static_cast<uint32_t>(lo % d);
        }

        carry = rem_u32_step(hi, carry, d);
        return rem_u32_step(lo, carry, d);
    } else {
        using Half = decltype(get_hi(val));
        Half hi = get_hi(val);
        Half lo = get_lo(val);

        if (carry == 0 && hi == Half(0)) {
            return rem_u32_step_generic(lo, 0, d);
        }

        carry = rem_u32_step_generic(hi, carry, d);
        return rem_u32_step_generic(lo, carry, d);
    }
}

template <typename T>
inline uint32_t rem_u32(const T& val, uint32_t d) {
    return rem_u32_step_generic(val, 0, d);
}

inline uint32_t div_mod_u32_step(uint64_t val, uint32_t carry, uint32_t d, uint64_t& out_q) {
    uint32_t hi = static_cast<uint32_t>(val >> 32);
    uint32_t lo = static_cast<uint32_t>(val);

    uint64_t cur1 = (static_cast<uint64_t>(carry) << 32) | hi;
    uint32_t q_hi = static_cast<uint32_t>(cur1 / d);
    carry = static_cast<uint32_t>(cur1 % d);

    uint64_t cur2 = (static_cast<uint64_t>(carry) << 32) | lo;
    uint32_t q_lo = static_cast<uint32_t>(cur2 / d);
    carry = static_cast<uint32_t>(cur2 % d);

    out_q = (static_cast<uint64_t>(q_hi) << 32) | q_lo;
    return carry;
}

template <typename T>
inline uint32_t div_mod_u32_step_generic(const T& val, uint32_t carry, uint32_t d, T& out_q) {
    if constexpr (std::is_integral_v<T>) {
        uint64_t q64 = 0;
        carry = div_mod_u32_step(static_cast<uint64_t>(val), carry, d, q64);
        out_q = static_cast<T>(q64);
        return carry;
    } else if constexpr (std::is_same_v<T, ::bignum::u128::U128>) {
        uint64_t hi = get_hi(val);
        uint64_t lo = get_lo(val);
        uint64_t q_hi = 0, q_lo = 0;

        if (carry == 0 && hi == 0) {
            carry = div_mod_u32_step(lo, 0, d, q_lo);
            out_q = make_big<::bignum::u128::U128>(uint64_t(0), q_lo);
            return carry;
        }

        carry = div_mod_u32_step(hi, carry, d, q_hi);
        carry = div_mod_u32_step(lo, carry, d, q_lo);
        out_q = make_big<::bignum::u128::U128>(q_hi, q_lo);
        return carry;
    } else {
        using Half = decltype(get_hi(val));
        Half hi = get_hi(val);
        Half lo = get_lo(val);
        Half q_hi{}, q_lo{};

        if (carry == 0 && hi == Half(0)) {
            carry = div_mod_u32_step_generic(lo, 0, d, q_lo);
            out_q = make_big<T>(Half(0), q_lo);
            return carry;
        }

        carry = div_mod_u32_step_generic(hi, carry, d, q_hi);
        carry = div_mod_u32_step_generic(lo, carry, d, q_lo);
        out_q = make_big<T>(q_hi, q_lo);
        return carry;
    }
}

template <typename T>
inline T div_by_u32(const T& val, uint32_t d) {
    T q{};
    div_mod_u32_step_generic(val, 0, d, q);
    return q;
}

// --- Универсальные операции общего деления UBig ---

template <typename T>
inline constexpr std::pair<T, T> div_rem(const T& a, const T& b) {
    if constexpr (requires { { a / b } -> std::same_as<std::pair<T, T>>; }) {
        return a / b;
    } else if constexpr (requires { T::template divide<true, true>(a, b, (T*)nullptr); }) {
        T r;
        T q = T::template divide<true, true>(a, b, &r);
        return {q, r};
    } else {
        return { a / b, a % b };
    }
}

template <typename T>
inline constexpr T div_q(const T& a, const T& b) {
    return div_rem(a, b).first;
}

template <typename T>
inline constexpr T div_r(const T& a, const T& b) {
    return div_rem(a, b).second;
}

template <typename T>
inline T gcd_generic(T x, T y) {
    while (y != T(0)) {
        T r = div_r(x, y);
        x = y;
        y = r;
    }
    return x;
}

template <typename T>
inline void mult_mod(T& x, const T& y, const T& m) {
    using DoubleT = bignum::UBig<T>;
    DoubleT prod = DoubleT::mult_ext(x, y);
    x = (prod / m).second;
}

template <typename T>
inline void square_mod(T& x, const T& m) {
    using DoubleT = bignum::UBig<T>;
    DoubleT sq = DoubleT::square_ext(x);
    x = (sq / m).second;
}

template <typename T>
inline void int_power_mod(T& x, T exp, const T& m) {
    T base = x;
    x = T(1);
    while (exp > T(0)) {
        if ((exp & T(1)) == T(1)) mult_mod(x, base, m);
        exp >>= 1;
        square_mod(base, m);
    }
}

// --- Вспомогательная модульная арифметика ---

inline uint64_t mod_pow32(uint64_t base, uint32_t exp, uint32_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return res;
}

inline int legendre_symbol(uint32_t a, uint32_t p) {
    if (p == 2) return 1;
    uint32_t rem = a % p;
    if (rem == 0) return 0;
    uint64_t res = mod_pow32(rem, (p - 1) / 2, p);
    return (res == 1) ? 1 : -1;
}

inline uint32_t tonelli_shanks_u32(uint32_t n, uint32_t p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;
    if (legendre_symbol(n, p) != 1) return 0;
    if ((p & 3) == 3) return static_cast<uint32_t>(mod_pow32(n, (p + 1) / 4, p));

    uint32_t q = p - 1;
    uint32_t s = 0;
    while ((q & 1) == 0) { q >>= 1; s++; }

    uint32_t z = 2;
    while (legendre_symbol(z, p) != -1) z++;

    uint64_t c = mod_pow32(z, q, p);
    uint64_t r = mod_pow32(n, (q + 1) / 2, p);
    uint64_t t = mod_pow32(n, q, p);
    uint32_t m = s;

    while (t != 1) {
        uint32_t i = 1;
        uint64_t t2 = (t * t) % p;
        while (t2 != 1) {
            t2 = (t2 * t2) % p;
            i++;
        }
        uint64_t b = c;
        for (uint32_t j = 0; j < m - i - 1; ++j) b = (b * b) % p;
        r = (r * b) % p;
        c = (b * b) % p;
        t = (t * c) % p;
        m = i;
    }
    return static_cast<uint32_t>(r);
}

inline uint32_t mod_inv32(uint32_t a, uint32_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    int64_t a0 = a % m;
    if (m == 1) return 0;
    while (a0 > 1) {
        q = a0 / m0;
        t = m0;
        m0 = a0 % m0;
        a0 = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += m;
    return static_cast<uint32_t>(x1);
}

// --- Структуры данных ---

struct FactorBasePrime {
    uint32_t p;
    uint32_t r;     // r^2 = N mod p
    uint8_t log_p;
};

template <typename T>
struct SmoothRelation {
    T u;                              // (Ax + B) mod N
    std::vector<uint32_t> prime_idx; // Индексы простых из FB
    std::vector<uint16_t> prime_pow; // Степени
    bool is_negative = false;
    uint64_t large_prime = 1;         // Большой делитель L (1, если чисто гладкое)
};

// --- Основной класс Квадратичного Решета (SIQS) ---

template <typename T>
class MPQSSolver {
    T m_N;
    std::vector<FactorBasePrime> m_FB;
    int m_M = 32768;
    int m_num_relations_needed = 0;
    uint64_t m_L_max = 0;

    std::mutex m_rel_mutex;
    std::vector<SmoothRelation<T>> m_relations;
    std::unordered_set<T, BigHash<T>> m_seen_u;
    std::unordered_map<uint64_t, SmoothRelation<T>> m_partials;
    size_t m_combined_count = 0;
    std::atomic<bool> m_done{false};

public:
    explicit MPQSSolver(const T& n) : m_N(n) {}

    std::optional<T> solve(std::function<void(const std::string&)> status_cb = nullptr) {
        m_seen_u.clear();
        uint32_t bits = m_N.bit_width();
        if (bits > 320) {
            return std::nullopt;
        }

        // Адаптивный размер факторной базы
        size_t fb_size = 500;
        if (bits > 120) fb_size = 900;
        if (bits > 150) fb_size = 1600;
        if (bits > 180) fb_size = 3000;
        if (bits > 210) fb_size = 5000;
        if (bits > 240) fb_size = 8500;
        if (bits > 270) fb_size = 14000;

        m_M = (bits > 200) ? 65536 : 32768;
        m_num_relations_needed = static_cast<int>(fb_size) + 65;

        // 1. Построение факторной базы
        if (status_cb) status_cb("SIQS: Построение факторной базы...");
        if (!build_factor_base(fb_size)) {
            return std::nullopt;
        }

        m_relations.clear();
        m_partials.clear();
        m_combined_count = 0;
        m_done.store(false);

        // 2. Параллельное просеивание полиномов
        if (status_cb) status_cb("SIQS: Просеивание полиномов (SIQS Gray Code)...");
        unsigned num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;

        std::vector<std::thread> workers;
        workers.reserve(num_threads);

        std::atomic<uint32_t> batch_seed{0};

        for (unsigned t = 0; t < num_threads; ++t) {
            workers.emplace_back([this, &batch_seed, t, status_cb]() {
                sieve_worker(batch_seed, t, status_cb);
            });
        }

        for (auto& w : workers) {
            if (w.joinable()) w.join();
        }

        if (::u128::Globals::LoadStop() || m_relations.size() < static_cast<size_t>(m_FB.size() + 1)) {
            return std::nullopt;
        }

        // 3. Линейная алгебра над GF(2)
        if (status_cb) status_cb("SIQS: Решение СЛАУ (метод Гаусса над GF(2))...");
        return solve_linear_algebra(status_cb);
    }

private:
    bool build_factor_base(size_t target_size) {
        size_t prime_bound = std::max<size_t>(100000, target_size * 25);
        auto all_primes = ::u128::utils::primes(prime_bound);
        m_FB.clear();
        m_FB.reserve(target_size);

        m_FB.push_back({2, to_u32_low(m_N) & 1, 1});

        for (unsigned p : all_primes) {
            if (p == 2) continue;
            if (::u128::Globals::LoadStop()) return false;

            uint32_t n_mod_p = rem_u32(m_N, p);
            if (n_mod_p == 0) {
                return false;
            }

            if (legendre_symbol(n_mod_p, p) == 1) {
                uint32_t r = tonelli_shanks_u32(n_mod_p, p);
                uint8_t lp = static_cast<uint8_t>(std::round(std::log2(p)));
                m_FB.push_back({p, r, lp});
                if (m_FB.size() >= target_size) break;
            }
        }

        if (m_FB.size() < target_size) return false;

        uint64_t max_fb_p = m_FB.back().p;
        m_L_max = std::min<uint64_t>(max_fb_p * max_fb_p, max_fb_p * 256ULL);

        return true;
    }

    void sieve_worker(std::atomic<uint32_t>& batch_seed, unsigned thread_id,
                      std::function<void(const std::string&)> status_cb) {
        std::vector<uint8_t> sieve(2 * m_M, 0);
        size_t fb_len = m_FB.size();

        // Пропуск малых простых (<= 37) при просеивании
        size_t sieve_start_idx = 0;
        while (sieve_start_idx < fb_len && m_FB[sieve_start_idx].p <= 37) {
            sieve_start_idx++;
        }

        // Оптимальный битовый размер A: A ~ sqrt(2N) / M
        uint32_t half_bits = m_N.bit_width() / 2;
        uint32_t m_log2 = static_cast<uint32_t>(std::round(std::log2(m_M)));
        uint32_t target_a_bits = (half_bits > m_log2) ? (half_bits - m_log2) : 24;

        // Запас для порога: log2(L_max) + 24 бита на пропущенные p <= 37 и вариацию Q(x)
        uint8_t l_max_bits = (m_L_max > 0) ? static_cast<uint8_t>(std::round(std::log2(m_L_max))) : 0;
        uint8_t allowed_slack = l_max_bits + 24;

        // Предвыделенные структуры для Gray Code
        std::vector<uint32_t> a_inv(fb_len, 0);
        std::vector<uint32_t> root1(fb_len, 0);
        std::vector<uint32_t> root2(fb_len, 0);
        std::vector<bool> is_in_a(fb_len, false);
        std::vector<std::vector<uint32_t>> deltas;

        std::vector<uint32_t> p_indices;
        std::vector<uint16_t> p_powers;
        p_indices.reserve(64);
        p_powers.reserve(64);

        std::mt19937_64 rng(thread_id * 10007 + 1);

        while (!m_done.load(std::memory_order_relaxed) && !::u128::Globals::LoadStop()) {
            uint32_t b_id = batch_seed.fetch_add(1);

            // Подбор числа простых k для Carrier Polynomials
            size_t k_primes = 4;
            if (target_a_bits > 45) k_primes = 5;
            if (target_a_bits > 65) k_primes = 6;
            if (target_a_bits > 85) k_primes = 7;
            if (target_a_bits > 110) k_primes = 8;

            // Выбираем простые ближе к оптимуму: target_a_bits / k_primes
            uint32_t prime_target_bits = target_a_bits / k_primes;
            size_t p_min_idx = sieve_start_idx;
            while (p_min_idx < fb_len && m_FB[p_min_idx].log_p < prime_target_bits - 1) {
                p_min_idx++;
            }
            size_t p_max_idx = p_min_idx;
            while (p_max_idx < fb_len && m_FB[p_max_idx].log_p <= prime_target_bits + 2) {
                p_max_idx++;
            }
            if (p_max_idx <= p_min_idx) {
                p_min_idx = fb_len / 4;
                p_max_idx = fb_len;
            }

            std::vector<size_t> a_indices;
            a_indices.reserve(k_primes);
            T A = T(1);

            std::fill(is_in_a.begin(), is_in_a.end(), false);
            size_t attempts = 0;
            while (a_indices.size() < k_primes && attempts < 200) {
                attempts++;
                size_t p_idx = p_min_idx + (rng() % (p_max_idx - p_min_idx));
                if (!is_in_a[p_idx]) {
                    is_in_a[p_idx] = true;
                    a_indices.push_back(p_idx);
                    A = A * T(m_FB[p_idx].p);
                }
            }
            if (a_indices.size() < k_primes || A == T(0)) continue;

            // Компоненты B_i через CRT
            std::vector<T> B_components(k_primes);
            for (size_t i = 0; i < k_primes; ++i) {
                uint32_t q = m_FB[a_indices[i]].p;
                uint32_t r = m_FB[a_indices[i]].r;
                T M_q = div_by_u32(A, q);
                uint32_t inv_q = mod_inv32(rem_u32(M_q, q), q);
                uint32_t gamma = static_cast<uint32_t>((1ULL * r * inv_q) % q);
                if (gamma > (q >> 1)) {
                    gamma = q - gamma; // Минимизируем размер компоненты
                }
                B_components[i] = M_q * T(gamma);
            }

            // Предвычисление инверсий A^-1 mod p и шагов корней deltas[j][p]
            deltas.assign(k_primes, std::vector<uint32_t>(fb_len, 0));
            for (size_t i = sieve_start_idx; i < fb_len; ++i) {
                if (is_in_a[i]) continue;
                uint32_t p = m_FB[i].p;
                uint32_t a_mod = rem_u32(A, p);
                if (a_mod == 0) continue;
                uint32_t inv = mod_inv32(a_mod, p);
                a_inv[i] = inv;

                for (size_t j = 1; j < k_primes; ++j) {
                    uint32_t b_mod = rem_u32(B_components[j], p);
                    deltas[j][i] = static_cast<uint32_t>((2ULL * b_mod * inv) % p);
                }
            }

            // Вектор текущих знаков компонент (+1 / -1)
            std::vector<int8_t> signs(k_primes, 1);

            // Начальные корни для базового B (все знаки +1)
            for (size_t i = sieve_start_idx; i < fb_len; ++i) {
                if (is_in_a[i]) continue;
                uint32_t p = m_FB[i].p;
                uint32_t r = m_FB[i].r;
                uint32_t inv = a_inv[i];

                uint32_t b_mod = 0;
                for (size_t k = 0; k < k_primes; ++k) {
                    b_mod = (b_mod + rem_u32(B_components[k], p)) % p;
                }

                root1[i] = static_cast<uint32_t>((1ULL * (p + r - b_mod) * inv) % p);
                root2[i] = static_cast<uint32_t>((1ULL * (2ULL * p - r - b_mod) * inv) % p);
            }

            size_t num_gray_polys = 1ULL << (k_primes - 1);

            // Перебор 2^(k-1) полиномов по коду Грея
            for (size_t g = 0; g < num_gray_polys; ++g) {
                if (m_done.load(std::memory_order_relaxed) || ::u128::Globals::LoadStop()) return;

                if (g > 0) {
                    uint32_t prev_gray = (g - 1) ^ ((g - 1) >> 1);
                    uint32_t cur_gray = g ^ (g >> 1);
                    uint32_t diff = prev_gray ^ cur_gray;
                    uint32_t v = ctz32(diff) + 1; // Компонента, меняющая знак
                    bool going_minus = (cur_gray & diff) != 0;

                    signs[v] = going_minus ? -1 : 1;

                    // Обновление корней СТРОГО в Z_p (без делений и редукций mod A!)
                    if (going_minus) {
                        for (size_t i = sieve_start_idx; i < fb_len; ++i) {
                            if (is_in_a[i]) continue;
                            uint32_t p = m_FB[i].p;
                            uint32_t d = deltas[v][i];
                            root1[i] = (root1[i] + d >= p) ? (root1[i] + d - p) : (root1[i] + d);
                            root2[i] = (root2[i] + d >= p) ? (root2[i] + d - p) : (root2[i] + d);
                        }
                    } else {
                        for (size_t i = sieve_start_idx; i < fb_len; ++i) {
                            if (is_in_a[i]) continue;
                            uint32_t p = m_FB[i].p;
                            uint32_t d = deltas[v][i];
                            root1[i] = (root1[i] >= d) ? (root1[i] - d) : (root1[i] + p - d);
                            root2[i] = (root2[i] >= d) ? (root2[i] - d) : (root2[i] + p - d);
                        }
                    }
                }

                // Точное вычисление |B| и знака B в целых числах без андерфлоу
                T B_pos = T(0);
                T B_neg = T(0);
                for (size_t k = 0; k < k_primes; ++k) {
                    if (signs[k] > 0) B_pos = B_pos + B_components[k];
                    else B_neg = B_neg + B_components[k];
                }

                T B_abs;
                bool B_is_negative = false;
                if (B_pos >= B_neg) {
                    B_abs = B_pos - B_neg;
                } else {
                    B_abs = B_neg - B_pos;
                    B_is_negative = true;
                }

                std::fill(sieve.begin(), sieve.end(), 0);

                // Просеивание массива
                for (size_t i = sieve_start_idx; i < fb_len; ++i) {
                    if (is_in_a[i]) continue;
                    uint32_t p = m_FB[i].p;
                    uint8_t log_p = m_FB[i].log_p;

                    int start1 = (static_cast<int>(root1[i]) + m_M) % p;
                    int start2 = (static_cast<int>(root2[i]) + m_M) % p;

                    for (int j = start1; j < 2 * m_M; j += p) sieve[j] += log_p;
                    for (int j = start2; j < 2 * m_M; j += p) sieve[j] += log_p;
                }

                // Оценка максимального значения |Q(x)| на отрезке [-M, M]
                T max_val = A * T(static_cast<uint64_t>(m_M) * m_M) + B_abs * T(static_cast<uint64_t>(2 * m_M)) + div_q(m_N, A);
                uint8_t max_bits = static_cast<uint8_t>(max_val.bit_width());
                uint8_t threshold = (max_bits > allowed_slack) ? (max_bits - allowed_slack) : 1;

                // Поиск кандидатов
                for (int j = 0; j < 2 * m_M; ++j) {
                    if (sieve[j] >= threshold) {
                        int64_t x_val = j - m_M;

                        // Вычисление Y = Ax + B с учетом знака B
                        T Ax = A * T(static_cast<uint64_t>(std::abs(x_val)));
                        T Y_abs;
                        bool Y_is_negative = false;

                        if (x_val >= 0) {
                            if (!B_is_negative) {
                                Y_abs = Ax + B_abs;
                            } else {
                                if (Ax >= B_abs) Y_abs = Ax - B_abs;
                                else { Y_abs = B_abs - Ax; Y_is_negative = true; }
                            }
                        } else {
                            if (!B_is_negative) {
                                if (B_abs >= Ax) Y_abs = B_abs - Ax;
                                else { Y_abs = Ax - B_abs; Y_is_negative = true; }
                            } else {
                                Y_abs = Ax + B_abs;
                                Y_is_negative = true;
                            }
                        }

                        // Q(x) = |Y^2 - N| / A
                        T Y2 = Y_abs * Y_abs;
                        T rem;
                        bool is_neg = false;

                        if (Y2 >= m_N) {
                            rem = div_q(Y2 - m_N, A);
                            is_neg = false;
                        } else {
                            rem = div_q(m_N - Y2, A);
                            is_neg = true;
                        }

                        if (rem == T(0)) continue;

                        p_indices.clear();
                        p_powers.clear();

                        // Быстрое извлечение двоек
                        if ((to_u32_low(rem) & 1) == 0) {
                            uint16_t c2 = 0;
                            while ((to_u32_low(rem) & 1) == 0 && rem > T(0)) {
                                rem >>= 1;
                                c2++;
                            }
                            p_indices.push_back(0);
                            p_powers.push_back(c2);
                        }

                        bool fits_u64 = (rem <= T(std::numeric_limits<uint64_t>::max()));
                        uint64_t rem64 = fits_u64 ? to_u64_low(rem) : 0;

                        // Факторизация пробными делениями
                        for (size_t pi = 1; pi < fb_len; ++pi) {
                            if (rem == T(1)) break;
                            uint32_t p = m_FB[pi].p;

                            if (fits_u64) {
                                if (rem64 == 1) break;
                                uint64_t p_sq = static_cast<uint64_t>(p) * p;
                                if (rem64 < p_sq) {
                                    if (rem64 <= m_FB.back().p) {
                                        auto it = std::lower_bound(m_FB.begin(), m_FB.end(), static_cast<uint32_t>(rem64),
                                                                   [](const FactorBasePrime& fbp, uint32_t val) { return fbp.p < val; });
                                        if (it != m_FB.end() && it->p == static_cast<uint32_t>(rem64)) {
                                            p_indices.push_back(static_cast<uint32_t>(std::distance(m_FB.begin(), it)));
                                            p_powers.push_back(1);
                                            rem64 = 1;
                                        }
                                    }
                                    break;
                                }

                                if (rem64 % p == 0) {
                                    uint16_t count = 0;
                                    do {
                                        count++;
                                        rem64 /= p;
                                    } while (rem64 % p == 0);
                                    p_indices.push_back(static_cast<uint32_t>(pi));
                                    p_powers.push_back(count);
                                }
                            } else {
                                uint32_t r = rem_u32(rem, p);
                                if (r == 0) {
                                    uint16_t count = 0;
                                    do {
                                        count++;
                                        rem = div_by_u32(rem, p);
                                        if (rem == T(1)) break;
                                        r = rem_u32(rem, p);
                                    } while (r == 0);
                                    p_indices.push_back(static_cast<uint32_t>(pi));
                                    p_powers.push_back(count);

                                    if (rem <= T(std::numeric_limits<uint64_t>::max())) {
                                        fits_u64 = true;
                                        rem64 = to_u64_low(rem);
                                    }
                                }
                            }
                        }

                        if (fits_u64) {
                            rem = T(rem64);
                        }

                        bool is_smooth = (rem == T(1));
                        bool is_partial = (!is_smooth && rem <= T(m_L_max));

                        if (is_smooth || is_partial) {
                            // Вычисляем u(x) = (Ax + B) mod N
                            T u_val = div_r(Y_abs, m_N);
                            if (Y_is_negative && u_val != T(0)) {
                                u_val = m_N - u_val;
                            }

                            // Сборка степеней простых
                            std::vector<std::pair<uint32_t, uint16_t>> factors;
                            factors.reserve(p_indices.size() + a_indices.size());
                            for (size_t k = 0; k < p_indices.size(); ++k) {
                                factors.emplace_back(p_indices[k], p_powers[k]);
                            }
                            for (size_t a_idx : a_indices) {
                                factors.emplace_back(static_cast<uint32_t>(a_idx), 1);
                            }
                            std::sort(factors.begin(), factors.end(), [](const auto& a, const auto& b) {
                                return a.first < b.first;
                            });

                            std::vector<uint32_t> final_indices;
                            std::vector<uint16_t> final_powers;
                            for (const auto& f : factors) {
                                if (!final_indices.empty() && final_indices.back() == f.first) {
                                    final_powers.back() += f.second;
                                } else {
                                    final_indices.push_back(f.first);
                                    final_powers.push_back(f.second);
                                }
                            }

                            if (is_smooth) {
                                T norm_u = (u_val > (m_N >> 1)) ? (m_N - u_val) : u_val;

                                std::lock_guard<std::mutex> lock(m_rel_mutex);
                                if (!m_seen_u.insert(norm_u).second) {
                                    continue;
                                }

                                SmoothRelation<T> rel;
                                rel.is_negative = is_neg;
                                rel.u = u_val;
                                rel.prime_idx = std::move(final_indices);
                                rel.prime_pow = std::move(final_powers);
                                rel.large_prime = 1;
                                m_relations.push_back(std::move(rel));

                                if (status_cb && (m_relations.size() % 10 == 0 ||
                                                  m_relations.size() >= static_cast<size_t>(m_num_relations_needed))) {
                                    status_cb("SIQS: Сбор соотношений: " + std::to_string(m_relations.size()) +
                                              " / " + std::to_string(m_num_relations_needed));
                                }

                                if (static_cast<int>(m_relations.size()) >= m_num_relations_needed) {
                                    m_done.store(true, std::memory_order_relaxed);
                                    return;
                                }
                            } else { // is_partial
                                uint64_t L = to_u64_low(rem);
                                SmoothRelation<T> rel;
                                rel.is_negative = is_neg;
                                rel.u = u_val;
                                rel.prime_idx = std::move(final_indices);
                                rel.prime_pow = std::move(final_powers);
                                rel.large_prime = L;

                                SmoothRelation<T> match_rel;
                                bool found_match = false;
                                {
                                    std::lock_guard<std::mutex> lock(m_rel_mutex);
                                    auto it = m_partials.find(L);
                                    if (it == m_partials.end()) {
                                        m_partials.emplace(L, std::move(rel));
                                    } else {
                                        match_rel = std::move(it->second);
                                        m_partials.erase(it);
                                        found_match = true;
                                    }
                                }

                                if (found_match) {
                                    SmoothRelation<T> combined;
                                    combined.u = match_rel.u;
                                    mult_mod(combined.u, rel.u, m_N);

                                    T norm_comb_u = (combined.u > (m_N >> 1)) ? (m_N - combined.u) : combined.u;
                                    combined.is_negative = (match_rel.is_negative != rel.is_negative);
                                    combined.large_prime = L;

                                    std::vector<std::pair<uint32_t, uint16_t>> comb_factors;
                                    comb_factors.reserve(match_rel.prime_idx.size() + rel.prime_idx.size());
                                    for (size_t k = 0; k < match_rel.prime_idx.size(); ++k) {
                                        comb_factors.emplace_back(match_rel.prime_idx[k], match_rel.prime_pow[k]);
                                    }
                                    for (size_t k = 0; k < rel.prime_idx.size(); ++k) {
                                        comb_factors.emplace_back(rel.prime_idx[k], rel.prime_pow[k]);
                                    }
                                    std::sort(comb_factors.begin(), comb_factors.end(), [](const auto& a, const auto& b) {
                                        return a.first < b.first;
                                    });

                                    for (const auto& f : comb_factors) {
                                        if (!combined.prime_idx.empty() && combined.prime_idx.back() == f.first) {
                                            combined.prime_pow.back() += f.second;
                                        } else {
                                            combined.prime_idx.push_back(f.first);
                                            combined.prime_pow.push_back(f.second);
                                        }
                                    }

                                    std::lock_guard<std::mutex> lock(m_rel_mutex);
                                    if (!m_seen_u.insert(norm_comb_u).second) {
                                        continue;
                                    }

                                    m_relations.push_back(std::move(combined));
                                    m_combined_count++;

                                    if (status_cb && (m_relations.size() % 10 == 0 ||
                                                      m_relations.size() >= static_cast<size_t>(m_num_relations_needed))) {
                                        status_cb("SIQS: Сбор соотношений: " + std::to_string(m_relations.size()) +
                                                  " / " + std::to_string(m_num_relations_needed) +
                                                  " (пар: " + std::to_string(m_combined_count) + ")");
                                    }

                                    if (static_cast<int>(m_relations.size()) >= m_num_relations_needed) {
                                        m_done.store(true, std::memory_order_relaxed);
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (b_id % 5 == 0 && thread_id == 0 && status_cb) {
                status_cb("SIQS: Пачек полиномов: " + std::to_string(b_id) +
                          " | Соотношений: " + std::to_string(m_relations.size()) +
                          " / " + std::to_string(m_num_relations_needed));
            }
        }
    }

    std::optional<T> solve_linear_algebra(std::function<void(const std::string&)> status_cb = nullptr) {
        size_t num_primes = m_FB.size();
        size_t num_rels = m_relations.size();
        if (num_rels <= num_primes) {
            if (status_cb) status_cb("SIQS: Ошибка: недостаточно соотношений (" +
                          std::to_string(num_rels) + " <= " + std::to_string(num_primes) + ")");
            return std::nullopt;
        }

        size_t total_cols = num_primes + 1;
        size_t words = (total_cols + 63) / 64;
        size_t hist_words = (num_rels + 63) / 64;

        struct MatrixRow {
            std::vector<uint64_t> mask;
            std::vector<uint64_t> history;
        };

        std::vector<MatrixRow> matrix(num_rels);

        for (size_t r = 0; r < num_rels; ++r) {
            matrix[r].mask.assign(words, 0);
            matrix[r].history.assign(hist_words, 0);
            matrix[r].history[r / 64] |= (1ULL << (r % 64));

            if (m_relations[r].is_negative) {
                matrix[r].mask[0] |= 1ULL;
            }

            for (size_t i = 0; i < m_relations[r].prime_idx.size(); ++i) {
                uint32_t col = m_relations[r].prime_idx[i] + 1;
                if (m_relations[r].prime_pow[i] & 1) {
                    matrix[r].mask[col / 64] ^= (1ULL << (col % 64));
                }
            }
        }

        // Метод Гаусса над GF(2)
        size_t row = 0;
        for (size_t col = 0; col < total_cols && row < num_rels; ++col) {
            size_t col_w = col / 64;
            uint64_t col_mask = 1ULL << (col % 64);

            size_t pivot = row;
            while (pivot < num_rels && !(matrix[pivot].mask[col_w] & col_mask)) {
                pivot++;
            }
            if (pivot >= num_rels) continue;

            if (row != pivot) {
                std::swap(matrix[row], matrix[pivot]);
            }

            const auto& pivot_mask = matrix[row].mask;
            const auto& pivot_hist = matrix[row].history;

            for (size_t r = 0; r < num_rels; ++r) {
                if (r != row && (matrix[r].mask[col_w] & col_mask)) {
                    // Начинаем XOR с col_w, так как предыдущие слова уже занулены
                    for (size_t w = col_w; w < words; ++w) {
                        matrix[r].mask[w] ^= pivot_mask[w];
                    }
                    for (size_t hw = 0; hw < hist_words; ++hw) {
                        matrix[r].history[hw] ^= pivot_hist[hw];
                    }
                }
            }
            row++;
        }

        int null_rows_count = 0;

        for (size_t r = 0; r < num_rels; ++r) {
            bool all_zero = true;
            for (size_t w = 0; w < words; ++w) {
                if (matrix[r].mask[w] != 0) {
                    all_zero = false;
                    break;
                }
            }

            if (all_zero) {
                size_t rels_in_combination = 0;
                for (size_t rel_idx = 0; rel_idx < num_rels; ++rel_idx) {
                    if (matrix[r].history[rel_idx / 64] & (1ULL << (rel_idx % 64))) {
                        rels_in_combination++;
                    }
                }

                if (rels_in_combination < 2) {
                    continue;
                }

                null_rows_count++;
                T X = T(1);
                T Y = T(1);
                std::vector<uint64_t> combined_powers(num_primes, 0);

                for (size_t rel_idx = 0; rel_idx < num_rels; ++rel_idx) {
                    if (matrix[r].history[rel_idx / 64] & (1ULL << (rel_idx % 64))) {
                        mult_mod(X, m_relations[rel_idx].u, m_N);

                        if (m_relations[rel_idx].large_prime > 1) {
                            T lp(m_relations[rel_idx].large_prime);
                            mult_mod(Y, lp, m_N);
                        }

                        for (size_t pi = 0; pi < m_relations[rel_idx].prime_idx.size(); ++pi) {
                            combined_powers[m_relations[rel_idx].prime_idx[pi]] +=
                                m_relations[rel_idx].prime_pow[pi];
                        }
                    }
                }

                for (size_t p_idx = 0; p_idx < num_primes; ++p_idx) {
                    uint64_t half_pow = combined_powers[p_idx] / 2;
                    if (half_pow > 0) {
                        T p_val(m_FB[p_idx].p);
                        int_power_mod(p_val, T(half_pow), m_N);
                        mult_mod(Y, p_val, m_N);
                    }
                }

                T X2 = X;
                square_mod(X2, m_N);
                T Y2 = Y;
                square_mod(Y2, m_N);

                if (X2 != Y2) {
                    continue;
                }

                if (X == Y) {
                    continue;
                }

                T diff = (X >= Y) ? (X - Y) : (Y - X);
                T sum = div_r(X + Y, m_N);

                if (sum == T(0)) {
                    continue;
                }

                T f1 = gcd_generic(diff, m_N);
                if (f1 > T(1) && f1 < m_N) {
                    if (status_cb) status_cb("SIQS: Успех! Найден множитель: " + f1.toString());
                    return f1;
                }

                T f2 = gcd_generic(sum, m_N);
                if (f2 > T(1) && f2 < m_N) {
                    if (status_cb) status_cb("SIQS: Успех! Найден множитель: " + f2.toString());
                    return f2;
                }
            }
        }

        if (status_cb) {
            status_cb("SIQS: Найдено зависимостей: " + std::to_string(null_rows_count) +
                      ", но все дали тривиальные делители.");
        }

        return std::nullopt;
    }
};

template <typename T>
inline std::optional<T> quadratic_sieve(const T& n, std::function<void(const std::string&)> status_cb = nullptr) {
    MPQSSolver<T> solver(n);
    return solver.solve(status_cb);
}

} // namespace bignum::siqs
