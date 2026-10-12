#pragma once

#include <cstdint>
#include <bit>

#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) || defined(__clang__)
#include <x86intrin.h>
#endif

namespace bignum::intrinsics {

inline uint8_t add_carry_u64(uint8_t carry_in, uint64_t a, uint64_t b, uint64_t& out_res) {
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
    return _addcarry_u64(carry_in, a, b, reinterpret_cast<unsigned __int64*>(&out_res));
#elif defined(__GNUC__) || defined(__clang__)
    return _addcarry_u64(carry_in, a, b, reinterpret_cast<unsigned long long*>(&out_res));
#else
    uint64_t s1 = a + carry_in;
    uint8_t c1 = (s1 < a) ? 1 : 0;
    out_res = s1 + b;
    uint8_t c2 = (out_res < s1) ? 1 : 0;
    return c1 | c2;
#endif
}

inline uint8_t sub_borrow_u64(uint8_t borrow_in, uint64_t a, uint64_t b, uint64_t& out_res) {
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
    return _subborrow_u64(borrow_in, a, b, reinterpret_cast<unsigned __int64*>(&out_res));
#elif defined(__GNUC__) || defined(__clang__)
    return _subborrow_u64(borrow_in, a, b, reinterpret_cast<unsigned long long*>(&out_res));
#else
    uint64_t d1 = a - borrow_in;
    uint8_t b1 = (a < borrow_in) ? 1 : 0;
    out_res = d1 - b;
    uint8_t b2 = (d1 < b) ? 1 : 0;
    return b1 | b2;
#endif
}

inline void mul64_to_128(uint64_t a, uint64_t b, uint64_t& out_hi, uint64_t& out_lo) {
#if defined(__SIZEOF_INT128__)
    unsigned __int128 prod = static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
    out_lo = static_cast<uint64_t>(prod);
    out_hi = static_cast<uint64_t>(prod >> 64);
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
    out_lo = _umul128(a, b, &out_hi);
#else
    uint32_t a_lo = static_cast<uint32_t>(a);
    uint32_t a_hi = static_cast<uint32_t>(a >> 32);
    uint32_t b_lo = static_cast<uint32_t>(b);
    uint32_t b_hi = static_cast<uint32_t>(b >> 32);

    uint64_t p0 = static_cast<uint64_t>(a_lo) * b_lo;
    uint64_t p1 = static_cast<uint64_t>(a_lo) * b_hi;
    uint64_t p2 = static_cast<uint64_t>(a_hi) * b_lo;
    uint64_t p3 = static_cast<uint64_t>(a_hi) * b_hi;

    uint64_t mid = p1 + (p0 >> 32);
    uint64_t mid_carry = (mid < p1) ? (1ULL << 32) : 0ULL;
    mid += p2;
    if (mid < p2) mid_carry += (1ULL << 32);

    out_lo = (mid << 32) | (p0 & 0xFFFFFFFFULL);
    out_hi = p3 + (mid >> 32) + mid_carry;
#endif
}

// Точное 128 x 128 -> 256 бит умножение
inline void mul128_ext(uint64_t a_hi, uint64_t a_lo,
                       uint64_t b_hi, uint64_t b_lo,
                       uint64_t& r0, uint64_t& r1, uint64_t& r2, uint64_t& r3) {
#if defined(__SIZEOF_INT128__)
    unsigned __int128 p0 = static_cast<unsigned __int128>(a_lo) * b_lo;
    unsigned __int128 p1 = static_cast<unsigned __int128>(a_lo) * b_hi;
    unsigned __int128 p2 = static_cast<unsigned __int128>(a_hi) * b_lo;
    unsigned __int128 p3 = static_cast<unsigned __int128>(a_hi) * b_hi;

    r0 = static_cast<uint64_t>(p0);

    unsigned __int128 mid = (p0 >> 64) + static_cast<uint64_t>(p1) + static_cast<uint64_t>(p2);
    r1 = static_cast<uint64_t>(mid);

    unsigned __int128 hi = (mid >> 64) + (p1 >> 64) + (p2 >> 64) + p3;
    r2 = static_cast<uint64_t>(hi);
    r3 = static_cast<uint64_t>(hi >> 64);
#else
    uint64_t l0, h0, l1, h1, l2, h2, l3, h3;
    mul64_to_128(a_lo, b_lo, h0, l0);
    mul64_to_128(a_lo, b_hi, h1, l1);
    mul64_to_128(a_hi, b_lo, h2, l2);
    mul64_to_128(a_hi, b_hi, h3, l3);

    r0 = l0;

    uint64_t mid1 = 0;
    uint8_t c1 = add_carry_u64(0, h0, l1, mid1);
    uint8_t c2 = add_carry_u64(0, mid1, l2, r1);
    uint64_t carry_to_r2 = static_cast<uint64_t>(c1) + c2;

    uint64_t mid2 = 0;
    uint8_t c3 = add_carry_u64(0, h1, h2, mid2);
    uint8_t c4 = add_carry_u64(0, mid2, l3, r2);
    uint8_t c5 = add_carry_u64(0, r2, carry_to_r2, r2);
    uint64_t carry_to_r3 = static_cast<uint64_t>(c3) + c4 + c5;

    r3 = h3 + carry_to_r3;
#endif
}

} // namespace bignum::intrinsics
