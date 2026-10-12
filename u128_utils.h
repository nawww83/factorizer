#pragma once

#include <atomic>
#include <vector>
#include <cstdint>

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