#pragma once

// ---------------------------------------------------------------------------
// siqs/common.hpp — общие типы решета и модулярная арифметика.
//
// Псевдонимы u32/u64/i64/Big, глобальный флаг остановки, 64-битная модулярка
// (умножение, степень, обратный элемент, простота, разложение на два
// множителя, корень Тонелли-Шенкса) и те же операции над Big.
//
// Внимание: системные заголовки подключаются до открытия namespace.
// ---------------------------------------------------------------------------

#include "bignum.hpp"
#include "cpu_topology.hpp"
#include "u128_utils.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <mutex>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(__SSE2__)
#  include <emmintrin.h>
#endif

namespace bignum::siqs_xor {

using u32 = uint32_t;
using u64 = uint64_t;
using i64 = int64_t;
using Big = UInt<512>;

inline bool stop_requested()
{
    return ::u128::Globals::LoadStop();
}

// --- 64-битная модулярная арифметика ---------------------------------------

inline u64 mul_mod64(u64 a, u64 b, u64 m)
{
#if defined(__SIZEOF_INT128__)
    return static_cast<u64>((static_cast<unsigned __int128>(a) * b) % m);
#else
    u64 hi, lo;
    intrinsics::mul64_to_128(a, b, hi, lo);
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
    u64 rem = 0;
    _udiv128(hi, lo, m, &rem);
    return rem;
#else
    u64 rem = 0;
    for (int bit = 127; bit >= 0; --bit)
    {
        const u64 bit_val = (bit >= 64) ? ((hi >> (bit - 64)) & 1ULL) : ((lo >> bit) & 1ULL);
        const u64 two = (rem << 1) | bit_val;
        if ((rem >> 63) || two >= m)
            two -= m;
        rem = two;
    }
    return rem;
#endif
#endif
}

inline u64 pow_mod64(u64 base, u64 exp, u64 m)
{
    if (m == 1)
        return 0;
    u64 res = 1;
    base %= m;
    while (exp > 0)
    {
        if (exp & 1)
            res = mul_mod64(res, base, m);
        base = mul_mod64(base, base, m);
        exp >>= 1;
    }
    return res;
}

inline u64 gcd64(u64 a, u64 b)
{
    while (b)
    {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

inline u64 inv_mod64(u64 a, u64 m)
{
    i64 t = 0, newt = 1;
    i64 r = static_cast<i64>(m), newr = static_cast<i64>(a % m);
    while (newr != 0)
    {
        i64 q = r / newr;
        i64 tmp_t = t - q * newt;
        t = newt;
        newt = tmp_t;
        i64 tmp_r = r - q * newr;
        r = newr;
        newr = tmp_r;
    }
    if (t < 0)
        t += static_cast<i64>(m);
    return static_cast<u64>(t);
}

// Детерминированный тест Миллера — Рабина для u64.
inline bool is_prime64(u64 n)
{
    if (n < 2)
        return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL,
                  23ULL, 29ULL, 31ULL, 37ULL})
    {
        if (n == p)
            return true;
        if (n % p == 0)
            return false;
    }
    u64 d = n - 1;
    const int r = static_cast<int>(std::countr_zero(d));
    d >>= r;
    for (u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL,
                  23ULL, 29ULL, 31ULL, 37ULL})
    {
        if (n <= a)
            break;
        u64 x = pow_mod64(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool composite = true;
        for (int i = 0; i < r - 1; ++i)
        {
            x = mul_mod64(x, x, n);
            if (x == n - 1)
            {
                composite = false;
                break;
            }
        }
        if (composite)
            return false;
    }
    return true;
}

// Метод Брента — Полларда для 64-битного кофактора.
inline bool brent_factor_two_primes(u64 n, u64 &f1, u64 &f2)
{
    if ((n & 1) == 0)
    {
        f1 = 2;
        f2 = n >> 1;
        return true;
    }
    for (u64 c : {1ULL, 3ULL, 5ULL})
    {
        u64 y = 2, m = 32, g = 1, r = 1, q = 1, ys = 2, x = 2;
        while (g == 1)
        {
            x = y;
            for (u64 i = 0; i < r; ++i)
                y = (mul_mod64(y, y, n) + c) % n;
            u64 k = 0;
            while (k < r && g == 1)
            {
                ys = y;
                const u64 lim = std::min(m, r - k);
                for (u64 i = 0; i < lim; ++i)
                {
                    y = (mul_mod64(y, y, n) + c) % n;
                    const u64 diff = (x > y) ? (x - y) : (y - x);
                    q = mul_mod64(q, diff, n);
                }
                g = gcd64(q, n);
                k += m;
            }
            r <<= 1;
            if (r > 1024)
                break;
        }
        if (g == n)
        {
            while (true)
            {
                ys = (mul_mod64(ys, ys, n) + c) % n;
                const u64 diff = (x > ys) ? (x - ys) : (ys - x);
                g = gcd64(diff, n);
                if (g > 1)
                    break;
            }
        }
        if (g > 1 && g < n)
        {
            f1 = g;
            f2 = n / g;
            return true;
        }
    }
    return false;
}

// Алгоритм Тонелли — Шэнкса: r^2 == n (mod p).
inline u32 tonelli_shanks(u64 n, u32 p)
{
    n %= p;
    if (n == 0)
        return 0;
    if (p == 2)
        return static_cast<u32>(n);
    if ((p & 3) == 3)
        return static_cast<u32>(pow_mod64(n, (p + 1) / 4, p));
    u32 q = p - 1, s = 0;
    while ((q & 1) == 0)
    {
        q >>= 1;
        ++s;
    }
    u32 z = 2;
    while (pow_mod64(z, (p - 1) / 2, p) != p - 1)
        ++z;
    u32 m = s;
    u64 c = pow_mod64(z, q, p);
    u64 t = pow_mod64(n, q, p);
    u64 r = pow_mod64(n, (q + 1) / 2, p);
    while (t != 1)
    {
        u32 i = 1;
        u64 temp = (t * t) % p;
        while (temp != 1)
        {
            temp = (temp * temp) % p;
            ++i;
        }
        const u64 b = pow_mod64(c, 1ULL << (m - i - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }
    return static_cast<u32>(r);
}

// --- Модулярные операции над 512-битным типом ------------------------------

inline Big pow_mod_big(Big base, u64 exp, const Big &m)
{
    BarrettReducer<512> br(m);
    Big res(1);
    if (base >= m)
        base = base % m;
    while (exp > 0)
    {
        if (exp & 1)
            res = br.mul(res, base);
        exp >>= 1;
        if (exp != 0)
            base = br.sqr(base);
    }
    return res;
}

inline Big inv_mod_big(Big a, const Big &m)
{
    a = a % m;
    Big u = a, v = m;
    Big x1(1), x2(0);
    while (!u.is_zero())
    {
        Big q, r;
        q = Big::divide(v, u, &r);
        v = u;
        u = r;
        Big qx1 = Big::mul_mod(q, x1, m);
        Big nx = (x2 >= qx1) ? (x2 - qx1) : (m - (qx1 - x2));
        x2 = x1;
        x1 = nx;
    }
    return x2;
}


} // namespace bignum::siqs_xor
