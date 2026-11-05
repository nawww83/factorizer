#pragma once

// ---------------------------------------------------------------------------
// siqs/config.hpp — параметры одного запуска решета.
//
// Всё, что считается один раз перед сбором соотношений: база простых,
// корни, границы больших простых, размер блока и пул простых для A.
// ---------------------------------------------------------------------------

#include "common.hpp"

namespace bignum::siqs_xor {

// --- Параметры и конфигурация ----------------------------------------------

struct SIQSConfig
{
    Big M;
    Big kM;
    std::vector<u32> prime_base;
    std::unordered_map<u64, u32> prime_to_idx;
    std::vector<u32> root1;
    std::vector<u32> root2;
    std::vector<uint8_t> log_primes;
    u32 num_primes = 0;
    u64 last_prime = 0;
    u64 min_2lp_comp = 0;
    u64 lp_bound = 0;
    u64 lp2_bound = 0;
    u32 B_size = 0;
    u32 half_B = 0;
    int s = 0;
    std::vector<u32> pool_indices;
    u32 small_prime_skip = 0;
    u32 split_idx = 0;
};

} // namespace bignum::siqs_xor
