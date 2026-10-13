#pragma once

#include <atomic>
#include <vector>
#include <cstdint>

#include "u128.hpp"

namespace u128
{

class Globals {
    static struct _gu128 {
        std::atomic<bool> is_stop = false;
    } global_u128;
public:
    explicit Globals() = default;
    static void SetStop(bool value) {
        global_u128.is_stop.store(value);
    }
    static bool LoadStop() {
        return global_u128.is_stop.load(std::memory_order::relaxed);
    }
};

namespace utils
{

// --- Универсальное извлечение младших бит из любой иерархии UBig ---

template<typename U>
inline constexpr uint64_t to_u64_low(const U &val)
{
    if constexpr (std::is_integral_v<U>) {
        return static_cast<uint64_t>(val);
    } else if constexpr (requires { val.low(); }) {
        return to_u64_low(val.low());
    } else {
        return 0;
    }
}

template<typename U>
inline constexpr uint32_t to_u32_low(const U &val)
{
    return static_cast<uint32_t>(to_u64_low(val));
}

// --- Извлечение половинок и сборка больших чисел ---

template<typename U>
inline constexpr auto get_hi(const U &val)
{
    if constexpr (requires { val.high(); }) {
        return val.high();
    } else if constexpr (requires { val.hi(); }) {
        return val.hi();
    } else if constexpr (std::is_same_v<U, U128>) {
        return to_u64_low(val >> 64);
    }
}

template<typename U>
inline constexpr auto get_lo(const U &val)
{
    if constexpr (requires { val.low(); }) {
        return val.low();
    } else if constexpr (requires { val.lo(); }) {
        return val.lo();
    } else if constexpr (std::is_same_v<U, U128>) {
        return to_u64_low(val);
    }
}

template<typename T, typename Half>
inline T make_big(const Half &hi, const Half &lo)
{
    if constexpr (std::is_same_v<T, U128>) {
        static const bool hi_first = []() {
            U128 test_obj{1ULL, 2ULL};
            return (to_u64_low(test_obj) == 2ULL);
        }();
        if (hi_first) {
            return U128{hi, lo};
        } else {
            return U128{lo, hi};
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

// --- Решето Эратосфена ---
inline std::vector<unsigned> primes(unsigned n) {
    if (n < 2) return {};

    std::vector<unsigned> ps;
    ps.reserve(n / 6);
    std::vector<uint8_t> is_prime_vec(n + 1, 1);

    ps.push_back(2);
    for (unsigned p = 3; p <= n; p += 2) {
        if (is_prime_vec[p]) {
            ps.push_back(p);
            if (1ULL * p * p <= n) {
                for (unsigned i = p * p; i <= n; i += 2 * p) {
                    is_prime_vec[i] = 0;
                }
            }
        }
    }
    return ps;
}


} // namespace utils

inline Globals::_gu128 Globals::global_u128;

} // namespace u128
