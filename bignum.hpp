#pragma once

// ---------------------------------------------------------------------------
// bignum.hpp — фиксированные беззнаковые целые на массиве 64-битных слов.
//
// Пришёл на смену рекурсивному UBig<T> из ubig.hpp: тип UInt<Bits> хранит
// число как плоский массив из Bits/64 слов, что даёт предсказуемую скорость
// и одинаковое поведение на любой разрядности (128 … 8192 бит).
//
// Умножение/сложение/вычитание используют аппаратные переносы из
// u128_intrinsics.hpp, деление — алгоритм D Кнута на 32-битных лимбах
// (работает и без __int128, то есть и под MSVC).
// ---------------------------------------------------------------------------

#include "u128_intrinsics.hpp"

#include <bit>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>

namespace bignum {

template <uint32_t Bits>
class UInt;

namespace detail {

// Деление «широкого» числа на «широкое» (алгоритм D Кнута) в системе
// счисления 2^32. u изменяется на месте. v[n-1] != 0, n >= 1.
// MaxLimbs — максимальное число 32-битных лимбов (m, n <= MaxLimbs).
template <int MaxLimbs>
inline void divmnu32(uint32_t *q, uint32_t *r, uint32_t *u,
                     const uint32_t *v, int m, int n) noexcept
{
    if (n == 1)
    {
        uint64_t k = 0;
        for (int j = m - 1; j >= 0; --j)
        {
            uint64_t cur = (k << 32) | u[j];
            q[j] = static_cast<uint32_t>(cur / v[0]);
            k = cur % v[0];
        }
        r[0] = static_cast<uint32_t>(k);
        return;
    }

    constexpr uint64_t B = 1ULL << 32;
    const int s = static_cast<int>(std::countl_zero(v[n - 1]));

    uint32_t vn[MaxLimbs];
    for (int i = n - 1; i > 0; --i)
        vn[i] = (v[i] << s) | (s ? (v[i - 1] >> (32 - s)) : 0);
    vn[0] = v[0] << s;

    uint32_t un[MaxLimbs + 1];
    un[m] = s ? (u[m - 1] >> (32 - s)) : 0;
    for (int i = m - 1; i > 0; --i)
        un[i] = (u[i] << s) | (s ? (u[i - 1] >> (32 - s)) : 0);
    un[0] = u[0] << s;

    for (int j = m - n; j >= 0; --j)
    {
        uint64_t num = (static_cast<uint64_t>(un[j + n]) << 32) | un[j + n - 1];
        uint64_t qhat = num / vn[n - 1];
        uint64_t rhat = num % vn[n - 1];
        while (qhat >= B || qhat * vn[n - 2] > ((rhat << 32) | un[j + n - 2]))
        {
            --qhat;
            rhat += vn[n - 1];
            if (rhat >= B)
                break;
        }

        int64_t k = 0;
        for (int i = 0; i < n; ++i)
        {
            uint64_t p = qhat * vn[i];
            int64_t t = static_cast<int64_t>(un[i + j]) - k
                        - static_cast<int64_t>(static_cast<uint32_t>(p));
            un[i + j] = static_cast<uint32_t>(t);
            k = static_cast<int64_t>(p >> 32) - (t >> 32);
        }
        int64_t t2 = static_cast<int64_t>(un[j + n]) - k;
        un[j + n] = static_cast<uint32_t>(t2);
        q[j] = static_cast<uint32_t>(qhat);

        if (t2 < 0)
        {
            q[j] -= 1;
            uint32_t carry = 0;
            for (int i = 0; i < n; ++i)
            {
                uint64_t ss = static_cast<uint64_t>(un[i + j]) + vn[i] + carry;
                un[i + j] = static_cast<uint32_t>(ss);
                carry = static_cast<uint32_t>(ss >> 32);
            }
            un[j + n] = static_cast<uint32_t>(un[j + n] + carry);
        }
    }

    for (int i = 0; i < n; ++i)
        r[i] = (un[i] >> s) | (s ? (un[i + 1] << (32 - s)) : 0);
}

} // namespace detail

template <uint32_t Bits>
class UInt
{
    static_assert(Bits % 64 == 0, "UInt<Bits>: Bits must be a multiple of 64");
    static_assert(Bits >= 64, "UInt<Bits>: Bits must be at least 64");

public:
    static constexpr uint32_t BITS = Bits;
    static constexpr uint32_t WORDS = Bits / 64;

    uint64_t w[WORDS];

    constexpr UInt() : w{} {}

    constexpr UInt(uint64_t v) : w{}
    {
        w[0] = v;
    }
    constexpr UInt(uint32_t v) : w{}
    {
        w[0] = v;
    }
    constexpr UInt(int v) : w{}
    {
        w[0] = static_cast<uint64_t>(static_cast<int64_t>(v));
    }

    template <uint32_t OtherBits>
    explicit UInt(const UInt<OtherBits> &o) : w{}
    {
        constexpr uint32_t n = (OtherBits < Bits) ? UInt<OtherBits>::WORDS : WORDS;
        for (uint32_t i = 0; i < n; ++i)
            w[i] = o.w[i];
    }

    UInt(const std::string &str) : w{}
    {
        *this = fromString(str);
    }

    // --- Состояние ---------------------------------------------------------

    [[nodiscard]] bool is_zero() const noexcept
    {
        for (uint32_t i = 0; i < WORDS; ++i)
            if (w[i] != 0)
                return false;
        return true;
    }

    [[nodiscard]] bool fits_u64() const noexcept
    {
        for (uint32_t i = 1; i < WORDS; ++i)
            if (w[i] != 0)
                return false;
        return true;
    }

    [[nodiscard]] uint64_t to_u64() const noexcept { return w[0]; }

    [[nodiscard]] uint64_t low_word() const noexcept { return w[0]; }

    [[nodiscard]] uint64_t &word(uint32_t i) noexcept { return w[i]; }
    [[nodiscard]] const uint64_t &word(uint32_t i) const noexcept { return w[i]; }

    [[nodiscard]] uint32_t bit_width() const noexcept
    {
        for (uint32_t i = WORDS; i-- > 0;)
            if (w[i] != 0)
                return i * 64 + static_cast<uint32_t>(std::bit_width(w[i]));
        return 0;
    }
    [[nodiscard]] uint32_t bit_length() const noexcept { return bit_width(); }

    [[nodiscard]] bool get_bit(uint32_t b) const noexcept
    {
        if (b >= Bits)
            return false;
        return ((w[b >> 6] >> (b & 63)) & 1ULL) != 0;
    }

    void set_bit(uint32_t b) noexcept
    {
        if (b < Bits)
            w[b >> 6] |= (1ULL << (b & 63));
    }

    // --- Сложение / вычитание (на месте) -----------------------------------

    void add(const UInt &o) noexcept
    {
        uint8_t carry = 0;
        for (uint32_t i = 0; i < WORDS; ++i)
            carry = intrinsics::add_carry_u64(carry, w[i], o.w[i], w[i]);
    }

    void sub(const UInt &o) noexcept
    {
        uint8_t borrow = 0;
        for (uint32_t i = 0; i < WORDS; ++i)
            borrow = intrinsics::sub_borrow_u64(borrow, w[i], o.w[i], w[i]);
    }

    void inc() noexcept { add(UInt(1)); }

    // --- Сдвиги (на месте) -------------------------------------------------

    void shl(uint32_t s) noexcept
    {
        if (s == 0)
            return;
        if (s >= Bits)
        {
            std::memset(w, 0, sizeof(w));
            return;
        }
        const uint32_t ws = s >> 6;
        const uint32_t bs = s & 63;
        if (ws)
        {
            for (uint32_t i = WORDS; i-- > ws;)
                w[i] = w[i - ws];
            for (uint32_t i = 0; i < ws; ++i)
                w[i] = 0;
        }
        if (bs)
        {
            for (uint32_t i = WORDS; i-- > 1;)
                w[i] = (w[i] << bs) | (w[i - 1] >> (64 - bs));
            w[0] <<= bs;
        }
    }

    void shr(uint32_t s) noexcept
    {
        if (s == 0)
            return;
        if (s >= Bits)
        {
            std::memset(w, 0, sizeof(w));
            return;
        }
        const uint32_t ws = s >> 6;
        const uint32_t bs = s & 63;
        if (ws)
        {
            for (uint32_t i = 0; i + ws < WORDS; ++i)
                w[i] = w[i + ws];
            for (uint32_t i = WORDS - ws; i < WORDS; ++i)
                w[i] = 0;
        }
        if (bs)
        {
            for (uint32_t i = 0; i + 1 < WORDS; ++i)
                w[i] = (w[i] >> bs) | (w[i + 1] << (64 - bs));
            w[WORDS - 1] >>= bs;
        }
    }

    void shl1() noexcept { shl(1); }
    void shr1() noexcept { shr(1); }

    // --- Деление на 64-битное число ---------------------------------------

    uint64_t divrem_u64(uint64_t d) noexcept
    {
        if (d == 0)
            return 0;
#if defined(__SIZEOF_INT128__)
        unsigned __int128 rem = 0;
        // Незначащие старшие слова пропускаем: числа в SIQS — это 3–6 слов, и
        // деление всех восьми слов на каждый простой базы стоило заметно.
        uint32_t i = WORDS;
        while (i > 0 && w[i - 1] == 0)
            --i;
        for (; i-- > 0;)
        {
            unsigned __int128 cur = (rem << 64) | w[i];
            w[i] = static_cast<uint64_t>(cur / d);
            rem = cur % d;
        }
        return static_cast<uint64_t>(rem);
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
        uint64_t rem = 0;
        uint32_t i = WORDS;
        while (i > 0 && w[i - 1] == 0)
            --i;
        for (; i-- > 0;)
        {
            uint64_t rem_out = 0;
            w[i] = _udiv128(rem, w[i], d, &rem_out);
            rem = rem_out;
        }
        return rem;
#else
        uint64_t rem = 0;
        const uint32_t top = bit_width();
        for (uint32_t b = top; b-- > 0;)
        {
            uint64_t two = (rem << 1) | static_cast<uint64_t>(get_bit(b));
            if ((rem >> 63) || two >= d)
                two -= d;
            rem = two;
        }
        return rem;
#endif
    }

    [[nodiscard]] uint64_t mod_u64(uint64_t d) const noexcept
    {
        UInt tmp = *this;
        return tmp.divrem_u64(d);
    }

    // --- Сравнения ---------------------------------------------------------

    friend bool operator==(const UInt &a, const UInt &b) noexcept
    {
        for (uint32_t i = 0; i < WORDS; ++i)
            if (a.w[i] != b.w[i])
                return false;
        return true;
    }
    friend bool operator!=(const UInt &a, const UInt &b) noexcept { return !(a == b); }

    friend bool operator<(const UInt &a, const UInt &b) noexcept
    {
        for (uint32_t i = WORDS; i-- > 0;)
        {
            if (a.w[i] < b.w[i])
                return true;
            if (a.w[i] > b.w[i])
                return false;
        }
        return false;
    }
    friend bool operator>(const UInt &a, const UInt &b) noexcept { return b < a; }
    friend bool operator<=(const UInt &a, const UInt &b) noexcept { return !(b < a); }
    friend bool operator>=(const UInt &a, const UInt &b) noexcept { return !(a < b); }

    // --- Арифметика --------------------------------------------------------

    friend UInt operator+(const UInt &a, const UInt &b) noexcept
    {
        UInt r = a;
        r.add(b);
        return r;
    }
    friend UInt operator-(const UInt &a, const UInt &b) noexcept
    {
        UInt r = a;
        r.sub(b);
        return r;
    }
    friend UInt operator&(const UInt &a, const UInt &b) noexcept
    {
        UInt r;
        for (uint32_t i = 0; i < WORDS; ++i)
            r.w[i] = a.w[i] & b.w[i];
        return r;
    }
    friend UInt operator|(const UInt &a, const UInt &b) noexcept
    {
        UInt r;
        for (uint32_t i = 0; i < WORDS; ++i)
            r.w[i] = a.w[i] | b.w[i];
        return r;
    }
    friend UInt operator^(const UInt &a, const UInt &b) noexcept
    {
        UInt r;
        for (uint32_t i = 0; i < WORDS; ++i)
            r.w[i] = a.w[i] ^ b.w[i];
        return r;
    }
    friend UInt operator~(const UInt &a) noexcept
    {
        UInt r;
        for (uint32_t i = 0; i < WORDS; ++i)
            r.w[i] = ~a.w[i];
        return r;
    }

    friend UInt operator<<(const UInt &a, uint32_t s) noexcept
    {
        UInt r = a;
        r.shl(s);
        return r;
    }
    friend UInt operator>>(const UInt &a, uint32_t s) noexcept
    {
        UInt r = a;
        r.shr(s);
        return r;
    }

    UInt &operator+=(const UInt &o) noexcept { add(o); return *this; }
    UInt &operator-=(const UInt &o) noexcept { sub(o); return *this; }
    UInt &operator&=(const UInt &o) noexcept
    {
        for (uint32_t i = 0; i < WORDS; ++i) w[i] &= o.w[i];
        return *this;
    }
    UInt &operator|=(const UInt &o) noexcept
    {
        for (uint32_t i = 0; i < WORDS; ++i) w[i] |= o.w[i];
        return *this;
    }
    UInt &operator^=(const UInt &o) noexcept
    {
        for (uint32_t i = 0; i < WORDS; ++i) w[i] ^= o.w[i];
        return *this;
    }
    UInt &operator<<=(uint32_t s) noexcept { shl(s); return *this; }
    UInt &operator>>=(uint32_t s) noexcept { shr(s); return *this; }

    UInt &operator++() noexcept { inc(); return *this; }
    UInt &operator--() noexcept { sub(UInt(1)); return *this; }
    UInt operator++(int) noexcept { UInt t = *this; inc(); return t; }
    UInt operator--(int) noexcept { UInt t = *this; sub(UInt(1)); return t; }

    // Умножение с усечением до Bits (как у обычного беззнакового типа).
    friend UInt operator*(const UInt &a, const UInt &b) noexcept
    {
        UInt r;
        for (uint32_t i = 0; i < WORDS; ++i)
        {
            if (a.w[i] == 0)
                continue;
            uint64_t carry = 0;
            for (uint32_t j = 0; i + j < WORDS; ++j)
            {
                uint64_t hi, lo;
                intrinsics::mul64_to_128(a.w[i], b.w[j], hi, lo);
                uint64_t t;
                uint8_t c1 = intrinsics::add_carry_u64(0, lo, r.w[i + j], t);
                uint8_t c2 = intrinsics::add_carry_u64(0, t, carry, r.w[i + j]);
                carry = hi + c1 + c2;
            }
        }
        return r;
    }

    UInt &operator*=(const UInt &o) noexcept { return *this = *this * o; }
    UInt &operator/=(const UInt &o) noexcept { return *this = divide(*this, o, nullptr); }
    UInt &operator%=(const UInt &o) noexcept
    {
        UInt r;
        divide(*this, o, &r);
        return *this = r;
    }

    friend UInt operator/(const UInt &a, const UInt &b) noexcept
    {
        return divide(a, b, nullptr);
    }
    friend UInt operator%(const UInt &a, const UInt &b) noexcept
    {
        UInt r;
        divide(a, b, &r);
        return r;
    }

    // --- Полное расширенное умножение (2 * Bits) ---------------------------

    static UInt<Bits * 2> mult_ext(const UInt &a, const UInt &b) noexcept
    {
        UInt<Bits * 2> r;
        constexpr uint32_t W2 = Bits * 2 / 64;
        for (uint32_t i = 0; i < WORDS; ++i)
        {
            if (a.w[i] == 0)
                continue;
            uint64_t carry = 0;
            for (uint32_t j = 0; j < WORDS; ++j)
            {
                uint64_t hi, lo;
                intrinsics::mul64_to_128(a.w[i], b.w[j], hi, lo);
                uint64_t t;
                uint8_t c1 = intrinsics::add_carry_u64(0, lo, r.w[i + j], t);
                uint8_t c2 = intrinsics::add_carry_u64(0, t, carry, r.w[i + j]);
                carry = hi + c1 + c2;
            }
            uint32_t k = i + WORDS;
            while (carry != 0 && k < W2)
            {
                carry = intrinsics::add_carry_u64(0, r.w[k], carry, r.w[k]);
                ++k;
            }
        }
        return r;
    }

    static UInt<Bits * 2> square_ext(const UInt &a) noexcept
    {
        return mult_ext(a, a);
    }

    template <uint32_t OutBits>
    [[nodiscard]] UInt<OutBits> truncate() const noexcept
    {
        UInt<OutBits> r;
        constexpr uint32_t n = (OutBits < Bits) ? UInt<OutBits>::WORDS : WORDS;
        for (uint32_t i = 0; i < n; ++i)
            r.w[i] = w[i];
        return r;
    }

    // --- Деление «широкое на широкое» -------------------------------------

    static UInt divide(const UInt &num, const UInt &den, UInt *rem_out = nullptr) noexcept
    {
        if (den.is_zero())
        {
            if (rem_out)
                *rem_out = UInt(0);
            return UInt(0);
        }
        if (num < den)
        {
            if (rem_out)
                *rem_out = num;
            return UInt(0);
        }
        if (num == den)
        {
            if (rem_out)
                *rem_out = UInt(0);
            return UInt(1);
        }
        if (num.fits_u64() && den.fits_u64())
        {
            const uint64_t n64 = num.w[0];
            const uint64_t d64 = den.w[0];
            if (rem_out)
                *rem_out = UInt(n64 % d64);
            return UInt(n64 / d64);
        }

        // Делитель помещается в одно слово — быстрый проход без Кнута.
        bool den_single_word = true;
        for (uint32_t i = 1; i < WORDS; ++i)
            if (den.w[i] != 0)
            {
                den_single_word = false;
                break;
            }
        if (den_single_word)
        {
            UInt q = num;
            uint64_t rem = q.divrem_u64(den.w[0]);
            if (rem_out)
                *rem_out = UInt(rem);
            return q;
        }

        // Общий случай: 32-битные лимбы, алгоритм D Кнута.
        constexpr int L = static_cast<int>(WORDS * 2);
        uint32_t u[L];
        uint32_t v[L];
        uint32_t ql[L];
        uint32_t rl[L];
        std::memset(u, 0, sizeof(u));
        std::memset(v, 0, sizeof(v));
        std::memset(ql, 0, sizeof(ql));
        std::memset(rl, 0, sizeof(rl));

        int m = 0, n = 0;
        for (uint32_t i = 0; i < WORDS; ++i)
        {
            u[2 * i] = static_cast<uint32_t>(num.w[i]);
            u[2 * i + 1] = static_cast<uint32_t>(num.w[i] >> 32);
            v[2 * i] = static_cast<uint32_t>(den.w[i]);
            v[2 * i + 1] = static_cast<uint32_t>(den.w[i] >> 32);
        }
        for (int i = static_cast<int>(WORDS * 2) - 1; i >= 0; --i)
            if (u[i] != 0)
            {
                m = i + 1;
                break;
            }
        for (int i = static_cast<int>(WORDS * 2) - 1; i >= 0; --i)
            if (v[i] != 0)
            {
                n = i + 1;
                break;
            }

        detail::divmnu32<L>(ql, rl, u, v, m, n);

        UInt q(0), r(0);
        for (int i = 0; i < n; ++i)
            r.w[i >> 1] |= static_cast<uint64_t>(rl[i]) << (32 * (i & 1));
        for (int i = 0; i < m; ++i)
            q.w[i >> 1] |= static_cast<uint64_t>(ql[i]) << (32 * (i & 1));

        if (rem_out)
            *rem_out = r;
        return q;
    }

    // (a * b) mod m — через полное произведение и точное деление.
    static UInt mul_mod(const UInt &a, const UInt &b, const UInt &m) noexcept
    {
        if (m.is_zero())
            return UInt(0);
        if (m == UInt(1))
            return UInt(0);
        const UInt ra = (a >= m) ? (a % m) : a;
        const UInt rb = (b >= m) ? (b % m) : b;
        const UInt<Bits * 2> prod = mult_ext(ra, rb);
        const UInt<Bits * 2> m2(m);
        const UInt<Bits * 2> rem = prod % m2;
        return rem.template truncate<Bits>();
    }

    // --- Строки ------------------------------------------------------------

    static UInt fromString(const std::string &str)
    {
        UInt result(0);
        size_t idx = 0;
        while (idx < str.size()
               && (str[idx] == ' ' || str[idx] == '\t' || str[idx] == '\n' || str[idx] == '\r'))
            ++idx;
        if (idx < str.size() && str[idx] == '+')
            ++idx;

        static constexpr uint64_t POW10[] = {
            1ULL, 10ULL, 100ULL, 1'000ULL, 10'000ULL, 100'000ULL,
            1'000'000ULL, 10'000'000ULL, 100'000'000ULL,
            1'000'000'000ULL, 10'000'000'000ULL, 100'000'000'000ULL,
            1'000'000'000'000ULL, 10'000'000'000'000ULL,
            100'000'000'000'000ULL, 1'000'000'000'000'000ULL,
            10'000'000'000'000'000ULL, 100'000'000'000'000'000ULL,
            1'000'000'000'000'000'000ULL};

        while (idx < str.size())
        {
            uint64_t chunk = 0;
            size_t count = 0;
            while (idx < str.size() && count < 18)
            {
                const char c = str[idx];
                if (c >= '0' && c <= '9')
                {
                    chunk = chunk * 10 + static_cast<uint64_t>(c - '0');
                    ++count;
                    ++idx;
                }
                else
                {
                    break;
                }
            }
            if (count == 0)
                break;
            result = result * UInt(POW10[count]) + UInt(chunk);
        }
        return result;
    }

    [[nodiscard]] std::string toString() const
    {
        if (is_zero())
            return "0";
        if (fits_u64())
            return std::to_string(w[0]);

        static constexpr uint64_t BASE18 = 1'000'000'000'000'000'000ULL;
        const UInt base(BASE18);
        std::string s;
        UInt cur = *this;
        while (!cur.is_zero())
        {
            UInt rem;
            cur = divide(cur, base, &rem);
            std::string chunk = std::to_string(rem.w[0]);
            if (!cur.is_zero() && chunk.size() < 18)
                chunk = std::string(18 - chunk.size(), '0') + chunk;
            s = chunk + s;
        }
        return s;
    }

    friend std::ostream &operator<<(std::ostream &os, const UInt &val)
    {
        return os << val.toString();
    }
};

// ---------------------------------------------------------------------------
// Алиасы разрядностей, используемых приложением
// ---------------------------------------------------------------------------

using U128 = UInt<128>;
using U256 = UInt<256>;
using U512 = UInt<512>;
using U1024 = UInt<1024>;
using U2048 = UInt<2048>;
using U4096 = UInt<4096>;

// ---------------------------------------------------------------------------
// Редуктор Барретта: быстрое умножение по модулю без деления
// ---------------------------------------------------------------------------

template <uint32_t Bits>
class BarrettReducer
{
    using T = UInt<Bits>;
    using DoubleT = UInt<Bits * 2>;

    T m_mod;
    DoubleT m_mu;
    uint32_t m_k = 0;

public:
    BarrettReducer() = default;

    explicit BarrettReducer(const T &mod) : m_mod(mod)
    {
        if (mod <= T(3))
            return;
        m_k = mod.bit_width();

        DoubleT half_pow(1);
        half_pow.shl(2 * m_k - 1);
        DoubleT den(mod);
        DoubleT r0;
        DoubleT q0 = DoubleT::divide(half_pow, den, &r0);

        DoubleT r_double = r0;
        r_double.shl(1);
        if (r_double >= den)
        {
            m_mu = q0;
            m_mu.shl(1);
            m_mu += DoubleT(1);
        }
        else
        {
            m_mu = q0;
            m_mu.shl(1);
        }

        // Проверка корректности редукции: (m-2)(m-3) mod m == 6, (m-1)^2 == 1.
        if (mod > T(6))
        {
            const T a = mod - T(2);
            const T b = mod - T(3);
            if (mul(a, b) != T(6) || sqr(mod - T(1)) != T(1))
                std::abort();
        }
    }

    [[nodiscard]] const T &modulus() const noexcept { return m_mod; }

    [[nodiscard]] inline T mul(const T &x, const T &y) const noexcept
    {
        return reduce(T::mult_ext(x, y));
    }

    [[nodiscard]] inline T sqr(const T &x) const noexcept
    {
        return reduce(T::square_ext(x));
    }

    [[nodiscard]] inline T reduce(const DoubleT &z) const noexcept
    {
        if (z < DoubleT(m_mod))
            return z.template truncate<Bits>();

        DoubleT q1 = z;
        q1.shr(m_k + 1);

        DoubleT q2 = q1 * m_mu;
        q2.shr(m_k - 1);

        DoubleT prod = q2 * DoubleT(m_mod);
        DoubleT r = (z >= prod) ? (z - prod) : DoubleT(0);

        while (r >= DoubleT(m_mod))
            r -= DoubleT(m_mod);

        return r.template truncate<Bits>();
    }
};

// ---------------------------------------------------------------------------
// Операции, общие для решета и ECM (бывший namespace bignum::siqs)
// ---------------------------------------------------------------------------

namespace ops {

template <uint32_t Bits>
inline UInt<Bits> div_q(const UInt<Bits> &a, const UInt<Bits> &b) noexcept
{
    return UInt<Bits>::divide(a, b, nullptr);
}

template <uint32_t Bits>
inline UInt<Bits> div_r(const UInt<Bits> &a, const UInt<Bits> &b) noexcept
{
    UInt<Bits> r;
    UInt<Bits>::divide(a, b, &r);
    return r;
}

template <uint32_t Bits>
inline std::pair<UInt<Bits>, UInt<Bits>> div_rem(const UInt<Bits> &a, const UInt<Bits> &b) noexcept
{
    UInt<Bits> r;
    UInt<Bits> q = UInt<Bits>::divide(a, b, &r);
    return {q, r};
}

template <uint32_t Bits>
inline UInt<Bits> gcd(UInt<Bits> a, UInt<Bits> b) noexcept
{
    while (!b.is_zero())
    {
        UInt<Bits> r = a % b;
        a = b;
        b = r;
    }
    return a;
}

template <uint32_t Bits>
inline uint32_t rem_u32(const UInt<Bits> &n, uint32_t d) noexcept
{
    return static_cast<uint32_t>(n.mod_u64(d));
}

template <uint32_t Bits>
inline UInt<Bits> div_by_u32(const UInt<Bits> &n, uint32_t d) noexcept
{
    UInt<Bits> q = n;
    q.divrem_u64(d);
    return q;
}

template <uint32_t Bits>
inline void mult_mod(UInt<Bits> &x, const UInt<Bits> &y, const UInt<Bits> &m) noexcept
{
    x = UInt<Bits>::mul_mod(x, y, m);
}

template <uint32_t Bits>
inline void square_mod(UInt<Bits> &x, const UInt<Bits> &m) noexcept
{
    x = UInt<Bits>::mul_mod(x, x, m);
}

template <uint32_t Bits>
inline UInt<Bits> pow_mod(UInt<Bits> base, UInt<Bits> exp, const UInt<Bits> &m) noexcept
{
    BarrettReducer<Bits> br(m);
    UInt<Bits> res(1);
    if (base >= m)
        base = base % m;
    while (!exp.is_zero())
    {
        if (exp.get_bit(0))
            res = br.mul(res, base);
        exp.shr(1);
        if (!exp.is_zero())
            base = br.sqr(base);
    }
    return res;
}

template <uint32_t Bits>
inline void int_power_mod(UInt<Bits> &x, UInt<Bits> exp, const UInt<Bits> &m) noexcept
{
    x = pow_mod(x, exp, m);
}

template <uint32_t Bits>
inline UInt<Bits> isqrt(const UInt<Bits> &n) noexcept
{
    if (n.is_zero())
        return UInt<Bits>(0);
    UInt<Bits> res(0);
    const int top = static_cast<int>((n.bit_width() - 1) >> 1);
    for (int i = top; i >= 0; --i)
    {
        UInt<Bits> cand = res;
        cand.set_bit(static_cast<uint32_t>(i));
        const UInt<Bits * 2> sq = UInt<Bits>::square_ext(cand);
        if (sq <= UInt<Bits * 2>(n))
            res = cand;
    }
    return res;
}

} // namespace ops

} // namespace bignum

// Глобальные алиасы — чтобы код приложения мог писать U512 без префикса.
using U128 = bignum::U128;
using U256 = bignum::U256;
using U512 = bignum::U512;
using U1024 = bignum::U1024;
using U2048 = bignum::U2048;
using U4096 = bignum::U4096;
