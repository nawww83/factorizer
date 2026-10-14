#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <bit>

namespace bignum::u128 {

class U128 {
    uint64_t m_lo = 0;
    uint64_t m_hi = 0;

public:
    using value_type = uint64_t;
    static constexpr uint32_t BITS = 128;

    constexpr U128() = default;
    constexpr U128(uint64_t lo) : m_lo(lo), m_hi(0) {}
    constexpr U128(uint32_t lo) : m_lo(lo), m_hi(0) {}
    constexpr U128(int lo) : m_lo(static_cast<uint64_t>(lo)), m_hi(0) {}
    // Порядок аргументов: (lo, hi), совпадающий с u128_utils.h
    constexpr U128(uint64_t lo, uint64_t hi) : m_lo(lo), m_hi(hi) {}

    explicit U128(const std::string& str) {
        *this = fromString(str);
    }

    constexpr uint64_t low()  const noexcept { return m_lo; }
    constexpr uint64_t high() const noexcept { return m_hi; }
    constexpr uint64_t lo()   const noexcept { return m_lo; }
    constexpr uint64_t hi()   const noexcept { return m_hi; }
    constexpr uint64_t to_u64() const noexcept { return m_lo; }

    // Проверка, помещается ли число в uint64_t
    [[nodiscard]] constexpr bool fits_u64() const noexcept {
        return m_hi == 0;
    }

    constexpr uint32_t bit_width() const {
        if (m_hi != 0) {
            return 128 - static_cast<uint32_t>(std::countl_zero(m_hi));
        }
        if (m_lo != 0) {
            return 64 - static_cast<uint32_t>(std::countl_zero(m_lo));
        }
        return 0;
    }

    // --- Сравнения ---

    friend constexpr bool operator==(const U128& a, const U128& b) {
        return a.m_lo == b.m_lo && a.m_hi == b.m_hi;
    }
    friend constexpr bool operator!=(const U128& a, const U128& b) {
        return !(a == b);
    }
    friend constexpr bool operator<(const U128& a, const U128& b) {
#if defined(__SIZEOF_INT128__)
        unsigned __int128 a_val = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        unsigned __int128 b_val = (static_cast<unsigned __int128>(b.m_hi) << 64) | b.m_lo;
        return a_val < b_val;
#else
        return (a.m_hi < b.m_hi) || (a.m_hi == b.m_hi && a.m_lo < b.m_lo);
#endif
    }
    friend constexpr bool operator>(const U128& a, const U128& b) {
        return b < a;
    }
    friend constexpr bool operator<=(const U128& a, const U128& b) {
        return !(b < a);
    }
    friend constexpr bool operator>=(const U128& a, const U128& b) {
        return !(a < b);
    }

    // --- Арифметика ---

    friend constexpr U128 operator+(const U128& a, const U128& b) {
#if defined(__SIZEOF_INT128__)
        unsigned __int128 a_val = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        unsigned __int128 b_val = (static_cast<unsigned __int128>(b.m_hi) << 64) | b.m_lo;
        unsigned __int128 sum = a_val + b_val;
        return U128(static_cast<uint64_t>(sum), static_cast<uint64_t>(sum >> 64));
#else
        uint64_t lo = 0, hi = 0;
        uint8_t c = intrinsics::add_carry_u64(0, a.m_lo, b.m_lo, lo);
        intrinsics::add_carry_u64(c, a.m_hi, b.m_hi, hi);
        return U128(lo, hi);
#endif
    }

    friend constexpr U128 operator-(const U128& a, const U128& b) {
#if defined(__SIZEOF_INT128__)
        unsigned __int128 a_val = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        unsigned __int128 b_val = (static_cast<unsigned __int128>(b.m_hi) << 64) | b.m_lo;
        unsigned __int128 diff = a_val - b_val;
        return U128(static_cast<uint64_t>(diff), static_cast<uint64_t>(diff >> 64));
#else
        uint64_t lo = 0, hi = 0;
        uint8_t borrow = intrinsics::sub_borrow_u64(0, a.m_lo, b.m_lo, lo);
        intrinsics::sub_borrow_u64(borrow, a.m_hi, b.m_hi, hi);
        return U128(lo, hi);
#endif
    }

    friend constexpr U128 operator*(const U128& a, const U128& b) {
#if defined(__SIZEOF_INT128__)
        unsigned __int128 a_val = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        unsigned __int128 b_val = (static_cast<unsigned __int128>(b.m_hi) << 64) | b.m_lo;
        unsigned __int128 prod = a_val * b_val;
        return U128(static_cast<uint64_t>(prod), static_cast<uint64_t>(prod >> 64));
#else
        uint64_t h0, l0;
        intrinsics::mul64_to_128(a.m_lo, b.m_lo, h0, l0);
        uint64_t hi = h0 + a.m_lo * b.m_hi + a.m_hi * b.m_lo;
        return U128(l0, hi);
#endif
    }

    // --- Побитовые операции ---

    friend constexpr U128 operator&(const U128& a, const U128& b) {
        return U128(a.m_lo & b.m_lo, a.m_hi & b.m_hi);
    }
    friend constexpr U128 operator|(const U128& a, const U128& b) {
        return U128(a.m_lo | b.m_lo, a.m_hi | b.m_hi);
    }
    friend constexpr U128 operator^(const U128& a, const U128& b) {
        return U128(a.m_lo ^ b.m_lo, a.m_hi ^ b.m_hi);
    }
    friend constexpr U128 operator~(const U128& a) {
        return U128(~a.m_lo, ~a.m_hi);
    }

    friend constexpr U128 operator<<(const U128& a, unsigned shift) {
        if (shift >= 128) return U128(0);
        if (shift == 0) return a;
#if defined(__SIZEOF_INT128__)
        unsigned __int128 v = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        v <<= shift;
        return U128(static_cast<uint64_t>(v), static_cast<uint64_t>(v >> 64));
#else
        if (shift >= 64) {
            return U128(0, a.m_lo << (shift - 64));
        }
        return U128(a.m_lo << shift, (a.m_hi << shift) | (a.m_lo >> (64 - shift)));
#endif
    }

    friend constexpr U128 operator>>(const U128& a, unsigned shift) {
        if (shift >= 128) return U128(0);
        if (shift == 0) return a;
#if defined(__SIZEOF_INT128__)
        unsigned __int128 v = (static_cast<unsigned __int128>(a.m_hi) << 64) | a.m_lo;
        v >>= shift;
        return U128(static_cast<uint64_t>(v), static_cast<uint64_t>(v >> 64));
#else
        if (shift >= 64) {
            return U128(a.m_hi >> (shift - 64), 0);
        }
        return U128((a.m_lo >> shift) | (a.m_hi << (64 - shift)), a.m_hi >> shift);
#endif
    }

    U128& operator++()    { return *this = *this + U128(1); }
    U128& operator--()    { return *this = *this - U128(1); }
    U128 operator++(int)  { U128 t = *this; ++(*this); return t; }
    U128 operator--(int)  { U128 t = *this; --(*this); return t; }

    U128& operator+=(const U128& o) { return *this = *this + o; }
    U128& operator-=(const U128& o) { return *this = *this - o; }
    U128& operator*=(const U128& o) { return *this = *this * o; }
    U128& operator&=(const U128& o) { return *this = *this & o; }
    U128& operator|=(const U128& o) { return *this = *this | o; }
    U128& operator^=(const U128& o) { return *this = *this ^ o; }
    U128& operator<<=(unsigned s)   { return *this = *this << s; }
    U128& operator>>=(unsigned s)   { return *this = *this >> s; }

    // --- Аппаратное деление 128 бит ---

    template <bool need_q = true, bool need_r = true>
    static U128 divide(const U128& num, const U128& den, U128* rem_out = nullptr) {
        if (den == U128(0)) {
            if (rem_out) *rem_out = U128(0);
            return U128(0);
        }
#if defined(__SIZEOF_INT128__)
        unsigned __int128 n = (static_cast<unsigned __int128>(num.m_hi) << 64) | num.m_lo;
        unsigned __int128 d = (static_cast<unsigned __int128>(den.m_hi) << 64) | den.m_lo;
        unsigned __int128 q = n / d;
        unsigned __int128 r = n % d;
        if constexpr (need_r) {
            if (rem_out) *rem_out = U128(static_cast<uint64_t>(r), static_cast<uint64_t>(r >> 64));
        }
        if constexpr (need_q) {
            return U128(static_cast<uint64_t>(q), static_cast<uint64_t>(q >> 64));
        } else {
            return U128(0);
        }
#else
        if (num < den) {
            if (rem_out) *rem_out = num;
            return U128(0);
        }
        if (num == den) {
            if (rem_out) *rem_out = U128(0);
            return U128(1);
        }

        uint32_t shift = num.bit_width() - den.bit_width();
        U128 current_den = den << shift;
        U128 rem = num;
        U128 q(0);

        for (int i = static_cast<int>(shift); i >= 0; --i) {
            if (rem >= current_den) {
                rem = rem - current_den;
                if constexpr (need_q) {
                    q = q | (U128(1) << i);
                }
            }
            current_den = current_den >> 1;
        }

        if constexpr (need_r) {
            if (rem_out) *rem_out = rem;
        }
        return q;
#endif
    }

    struct DivResult;

    friend DivResult operator/(const U128& a, const U128& b);
    friend U128 operator%(const U128& a, const U128& b);

    U128& operator/=(const U128& o);
    U128& operator%=(const U128& o);

    // --- Строки ---
    static U128 fromString(const std::string &str)
    {
        U128 result(0);
        size_t idx = 0;
        while (idx < str.size()
               && (str[idx] == ' ' || str[idx] == '\t' || str[idx] == '\n' || str[idx] == '\r')) {
            idx++;
        }
        if (idx < str.size() && str[idx] == '+') {
            idx++;
        }

        static constexpr uint64_t POW10[] = {1ULL,
                                             10ULL,
                                             100ULL,
                                             1'000ULL,
                                             10'000ULL,
                                             100'000ULL,
                                             1'000'000ULL,
                                             10'000'000ULL,
                                             100'000'000ULL,
                                             1'000'000'000ULL,
                                             10'000'000'000ULL,
                                             100'000'000'000ULL,
                                             1'000'000'000'000ULL,
                                             10'000'000'000'000ULL,
                                             100'000'000'000'000ULL,
                                             1'000'000'000'000'000ULL,
                                             10'000'000'000'000'000ULL,
                                             100'000'000'000'000'000ULL,
                                             1'000'000'000'000'000'000ULL};

        while (idx < str.size()) {
            uint64_t chunk = 0;
            size_t count = 0;
            while (idx < str.size() && count < 18) {
                char c = str[idx];
                if (c >= '0' && c <= '9') {
                    chunk = chunk * 10 + static_cast<uint64_t>(c - '0');
                    count++;
                    idx++;
                } else {
                    break;
                }
            }
            if (count == 0)
                break;
            result = result * U128(POW10[count]) + U128(chunk);
        }
        return result;
    }

    std::string toString() const
    {
        if (*this == U128(0))
            return "0";
        if (fits_u64()) {
            return std::to_string(to_u64());
        }

#if defined(__SIZEOF_INT128__)
        // Аппаратный путь через unsigned __int128:
        // Разбиваем 128 бит максимум на 3 блока по 10^18
        constexpr unsigned __int128 BASE18 = 1'000'000'000'000'000'000ULL;
        unsigned __int128 val = (static_cast<unsigned __int128>(m_hi) << 64) | m_lo;

        uint64_t r0 = static_cast<uint64_t>(val % BASE18);
        val /= BASE18;
        uint64_t r1 = static_cast<uint64_t>(val % BASE18);
        uint64_t r2 = static_cast<uint64_t>(val / BASE18);

        std::string s;
        if (r2 > 0) {
            s += std::to_string(r2);
            std::string s1 = std::to_string(r1);
            s += std::string(18 - s1.size(), '0') + s1;
            std::string s0 = std::to_string(r0);
            s += std::string(18 - s0.size(), '0') + s0;
        } else if (r1 > 0) {
            s += std::to_string(r1);
            std::string s0 = std::to_string(r0);
            s += std::string(18 - s0.size(), '0') + s0;
        } else {
            s += std::to_string(r0);
        }
        return s;
#else
        static constexpr uint64_t BASE18 = 1'000'000'000'000'000'000ULL;
        const U128 base_val(BASE18);
        std::string s;
        U128 cur = *this;

        while (cur > U128(0)) {
            if (cur.fits_u64()) {
                return std::to_string(cur.to_u64()) + s;
            }
            U128 rem(0);
            cur = divide<true, true>(cur, base_val, &rem);
            std::string chunk = std::to_string(rem.low());
            if (cur > U128(0) && chunk.size() < 18) {
                chunk = std::string(18 - chunk.size(), '0') + chunk;
            }
            s = chunk + s;
        }
        return s;
#endif
    }

    friend std::ostream &operator<<(std::ostream &os, const U128 &val)
    {
        return os << val.toString();
    }
};

struct U128::DivResult {
    U128 first;
    U128 second;

    constexpr operator U128() const { return first; }
};

inline U128::DivResult operator/(const U128& a, const U128& b) {
    U128 r(0);
    U128 q = U128::divide<true, true>(a, b, &r);
    return U128::DivResult{q, r};
}

inline U128 operator%(const U128& a, const U128& b) {
    U128 r(0);
    U128::divide<false, true>(a, b, &r);
    return r;
}

inline U128& U128::operator/=(const U128& o) {
    return *this = (*this / o).first;
}

inline U128& U128::operator%=(const U128& o) {
    return *this = *this % o;
}

} // namespace bignum::u128

namespace bignum {
using U128 = u128::U128;
}
using U128 = bignum::u128::U128;
