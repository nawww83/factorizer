#pragma once

#include <atomic>
#include <cstdint>
#include <vector>

#include "bignum.hpp"

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

// Младшие 64 бита любого целого типа (в т.ч. фиксированного UInt<Bits>).
template<typename U>
inline constexpr uint64_t to_u64_low(const U &val)
{
    if constexpr (std::is_integral_v<U>) {
        return static_cast<uint64_t>(val);
    } else if constexpr (requires { val.to_u64(); }) {
        return val.to_u64();
    } else {
        return 0;
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
