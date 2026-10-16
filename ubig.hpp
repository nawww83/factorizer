#pragma once

#include "u128.hpp"
#include "u128_intrinsics.hpp"
#include <bit>
#include <cstdint>
#include <string>
#include <iostream>
#include <type_traits>

namespace bignum {

template <typename Half>
class UBig {
    Half m_hi;
    Half m_lo;

public:
    using value_type = Half;
    using half_type  = Half;
    static constexpr uint32_t BITS = Half::BITS * 2;

    // Разрядность лимба для быстрого деления (алгоритм D Кнута): 64 бита там,
    // где есть 128-битные целые (GCC/Clang), иначе 32 бита (MSVC).
#if defined(__SIZEOF_INT128__)
    using limb_t = uint64_t;
    using wide_t = unsigned __int128; // двойной лимб
    using swide_t = __int128;         // знаковый двойной лимб
#else
    using limb_t = uint32_t;
    using wide_t = uint64_t;
    using swide_t = int64_t;
#endif
    static constexpr uint32_t LIMB_BITS = static_cast<uint32_t>(sizeof(limb_t) * 8);
    static constexpr uint32_t LIMBS = BITS / LIMB_BITS;
    static constexpr uint32_t HALF_LIMBS = Half::BITS / LIMB_BITS;
    static constexpr bool FLAT = (BITS % LIMB_BITS == 0) && (Half::BITS % LIMB_BITS == 0);

    constexpr UBig() : m_hi(Half(0)), m_lo(Half(0)) {}
    constexpr UBig(uint64_t val) : m_hi(Half(0)), m_lo(Half(val)) {}
    constexpr UBig(uint32_t val) : m_hi(Half(0)), m_lo(Half(val)) {}
    constexpr UBig(int val) : m_hi(Half(0)), m_lo(Half(static_cast<uint64_t>(val))) {}
    constexpr UBig(const Half& lo) : m_hi(Half(0)), m_lo(lo) {}
    constexpr UBig(const Half& hi, const Half& lo) : m_hi(hi), m_lo(lo) {}

    explicit UBig(const std::string& str) {
        *this = fromString(str);
    }

    constexpr operator Half() const noexcept { return m_lo; }

    constexpr const Half& high() const noexcept { return m_hi; }
    constexpr const Half& low()  const noexcept { return m_lo; }
    constexpr const Half& hi()   const noexcept { return m_hi; }
    constexpr const Half& lo()   const noexcept { return m_lo; }
    constexpr uint64_t to_u64() const noexcept { return m_lo.to_u64(); }

    [[nodiscard]] constexpr bool fits_u64() const noexcept {
        return (m_hi == Half(0)) && m_lo.fits_u64();
    }

    constexpr uint32_t bit_width() const {
        uint32_t w = m_hi.bit_width();
        if (w > 0) return Half::BITS + w;
        return m_lo.bit_width();
    }

    // --- Сравнения ---

    friend constexpr bool operator==(const UBig& a, const UBig& b) {
        return a.m_hi == b.m_hi && a.m_lo == b.m_lo;
    }
    friend constexpr bool operator!=(const UBig& a, const UBig& b) {
        return !(a == b);
    }
    friend constexpr bool operator<(const UBig& a, const UBig& b) {
        return (a.m_hi < b.m_hi) || (a.m_hi == b.m_hi && a.m_lo < b.m_lo);
    }
    friend constexpr bool operator>(const UBig& a, const UBig& b) {
        return b < a;
    }
    friend constexpr bool operator<=(const UBig& a, const UBig& b) {
        return !(b < a);
    }
    friend constexpr bool operator>=(const UBig& a, const UBig& b) {
        return !(a < b);
    }

    // --- Арифметика ---

    friend constexpr UBig operator+(const UBig& a, const UBig& b) {
        Half lo = a.m_lo + b.m_lo;
        Half carry = (lo < a.m_lo) ? Half(1) : Half(0);
        Half hi = a.m_hi + b.m_hi + carry;
        return UBig{hi, lo};
    }

    friend constexpr UBig operator-(const UBig& a, const UBig& b) {
        Half borrow = (a.m_lo < b.m_lo) ? Half(1) : Half(0);
        Half lo = a.m_lo - b.m_lo;
        Half hi = a.m_hi - b.m_hi - borrow;
        return UBig{hi, lo};
    }

    // --- Точное расширенное умножение (mult_ext) ---
    static UBig<Half> mult_ext(const Half& a, const Half& b) {
        if constexpr (std::is_same_v<Half, ::bignum::u128::U128>) {
            uint64_t r0, r1, r2, r3;
            ::bignum::intrinsics::mul128_ext(a.high(), a.low(), b.high(), b.low(), r0, r1, r2, r3);
            Half lo_part(r0, r1); // Для U128 конструктор принимает (lo, hi)
            Half hi_part(r2, r3);
            return UBig<Half>(hi_part, lo_part);
        } else {
            using Quarter = std::decay_t<decltype(a.low())>;
            using DoubleQuarter = UBig<Quarter>;

            DoubleQuarter z0 = DoubleQuarter::mult_ext(a.low(), b.low());
            DoubleQuarter z1 = DoubleQuarter::mult_ext(a.low(), b.high());
            DoubleQuarter z2 = DoubleQuarter::mult_ext(a.high(), b.low());
            DoubleQuarter z3 = DoubleQuarter::mult_ext(a.high(), b.high());

            Quarter r0 = z0.low();

            Quarter r1 = z0.high();
            Quarter c1 = Quarter(0);
            r1 += z1.low();
            if (r1 < z1.low()) c1 += Quarter(1);
            r1 += z2.low();
            if (r1 < z2.low()) c1 += Quarter(1);

            Quarter r2 = z1.high();
            Quarter c2 = Quarter(0);
            r2 += z2.high();
            if (r2 < z2.high()) c2 += Quarter(1);
            r2 += z3.low();
            if (r2 < z3.low()) c2 += Quarter(1);
            r2 += c1;
            if (r2 < c1) c2 += Quarter(1);

            Quarter r3 = z3.high() + c2;

            Half lo_res(r1, r0); // r1 — старшая четверть, r0 — младшая
            Half hi_res(r3, r2); // r3 — старшая четверть, r2 — младшая
            return UBig<Half>(hi_res, lo_res);
        }
    }

    static UBig<Half> square_ext(const Half& a) {
        return mult_ext(a, a);
    }

    friend constexpr UBig operator*(const UBig& a, const UBig& b) {
        UBig p0 = mult_ext(a.m_lo, b.m_lo);
        Half hi = p0.m_hi + a.m_lo * b.m_hi + a.m_hi * b.m_lo;
        return UBig{hi, p0.m_lo};
    }

    // --- Побитовые операции ---

    friend constexpr UBig operator&(const UBig& a, const UBig& b) {
        return UBig{a.m_hi & b.m_hi, a.m_lo & b.m_lo};
    }
    friend constexpr UBig operator|(const UBig& a, const UBig& b) {
        return UBig{a.m_hi | b.m_hi, a.m_lo | b.m_lo};
    }
    friend constexpr UBig operator^(const UBig& a, const UBig& b) {
        return UBig{a.m_hi ^ b.m_hi, a.m_lo ^ b.m_lo};
    }
    friend constexpr UBig operator~(const UBig& a) {
        return UBig{~a.m_hi, ~a.m_lo};
    }

    friend constexpr UBig operator<<(const UBig& a, unsigned shift) {
        if (shift >= BITS) return UBig{Half(0), Half(0)};
        if (shift == 0) return a;
        if (shift >= Half::BITS) {
            return UBig{a.m_lo << (shift - Half::BITS), Half(0)};
        }
        return UBig{(a.m_hi << shift) | (a.m_lo >> (Half::BITS - shift)), a.m_lo << shift};
    }

    friend constexpr UBig operator>>(const UBig& a, unsigned shift) {
        if (shift >= BITS) return UBig{Half(0), Half(0)};
        if (shift == 0) return a;
        if (shift >= Half::BITS) {
            return UBig{Half(0), a.m_hi >> (shift - Half::BITS)};
        }
        return UBig{a.m_hi >> shift, (a.m_lo >> shift) | (a.m_hi << (Half::BITS - shift))};
    }

    UBig& operator++()    { return *this = *this + UBig(1); }
    UBig& operator--()    { return *this = *this - UBig(1); }
    UBig operator++(int)  { UBig t = *this; ++(*this); return t; }
    UBig operator--(int)  { UBig t = *this; --(*this); return t; }

    UBig& operator+=(const UBig& o) { return *this = *this + o; }
    UBig& operator-=(const UBig& o) { return *this = *this - o; }
    UBig& operator*=(const UBig& o) { return *this = *this * o; }
    UBig& operator&=(const UBig& o) { return *this = *this & o; }
    UBig& operator|=(const UBig& o) { return *this = *this | o; }
    UBig& operator^=(const UBig& o) { return *this = *this ^ o; }
    UBig& operator<<=(unsigned s)   { return *this = *this << s; }
    UBig& operator>>=(unsigned s)   { return *this = *this >> s; }

    // --- Деление ---
    template<bool need_q = true, bool need_r = true>
    static UBig divide(const UBig &num, const UBig &den, UBig *rem_out = nullptr)
    {
        if (den == UBig(0)) {
            if (rem_out)
                *rem_out = UBig(0);
            return UBig(0);
        }
        if (num < den) {
            if (rem_out)
                *rem_out = num;
            return UBig(0);
        }
        if (num == den) {
            if (rem_out)
                *rem_out = UBig(0);
            return UBig(1);
        }

        // 1. Аппаратный путь для чисел до 64 бит
        if (num.fits_u64() && den.fits_u64()) {
            uint64_t n64 = num.to_u64();
            uint64_t d64 = den.to_u64();
            if constexpr (need_r) {
                if (rem_out)
                    *rem_out = UBig(n64 % d64);
            }
            return UBig(n64 / d64);
        }

        // 2. Рекурсивный спуск: если оба числа помещаются в Half,
        // сразу передаем в Half::divide (для U256 это аппаратный __int128 в U128)
        if (num.high() == Half(0) && den.high() == Half(0)) {
            Half r_half(0);
            Half q_half = Half::template divide<need_q, need_r>(num.low(),
                                                                den.low(),
                                                                need_r ? &r_half : nullptr);
            if constexpr (need_r) {
                if (rem_out)
                    *rem_out = UBig(r_half);
            }
            return UBig(q_half);
        }

        // 3. Быстрый путь: длинное деление (алгоритм D Кнута) "широкого" числа
        // на "половинку" (den.high() == 0). Это горячий путь mult_mod/square_mod.
        if constexpr (FLAT)
        {
            if (den.high() == Half(0))
            {
                limb_t u[LIMBS] = {};
                limb_t v[HALF_LIMBS] = {};
                limb_t q_limbs[LIMBS] = {};
                limb_t r_limbs[HALF_LIMBS] = {};

                UBig tmp = num;
                for (uint32_t i = 0; i < LIMBS; ++i)
                {
                    u[i] = static_cast<limb_t>(tmp.to_u64());
                    tmp >>= LIMB_BITS;
                }
                UBig dtmp = den;
                for (uint32_t i = 0; i < HALF_LIMBS; ++i)
                {
                    v[i] = static_cast<limb_t>(dtmp.to_u64());
                    dtmp >>= LIMB_BITS;
                }

                // divmnu требует v[n-1] != 0 — отбрасываем незначащие старшие лимбы делителя.
                int n = static_cast<int>(HALF_LIMBS);
                while (n > 1 && v[n - 1] == 0)
                    --n;

                divmnu(q_limbs, r_limbs, u, v, static_cast<int>(LIMBS), n);

                if constexpr (need_r)
                {
                    if (rem_out)
                    {
                        UBig R(0);
                        for (uint32_t i = 0; i < HALF_LIMBS; ++i)
                            R |= UBig(r_limbs[i]) << (LIMB_BITS * i);
                        *rem_out = R;
                    }
                }
                if constexpr (need_q)
                {
                    UBig Q(0);
                    for (uint32_t i = 0; i < LIMBS; ++i)
                        Q |= UBig(q_limbs[i]) << (LIMB_BITS * i);
                    return Q;
                }
                return UBig(0);
            }
        }

        // 4. Резервное побитовое деление (для делителя шире одной "половинки")
        uint32_t shift = num.bit_width() - den.bit_width();
        UBig current_den = den << shift;
        UBig rem = num;
        UBig q(0);

        for (int i = static_cast<int>(shift); i >= 0; --i) {
            if (rem >= current_den) {
                rem = rem - current_den;
                if constexpr (need_q) {
                    // Прямой битовый сдвиг: ставим бит без сдвига UBig(1) << i
                    if (i >= static_cast<int>(Half::BITS)) {
                        // Бит попадает в старшую половину
                        q.m_hi = q.m_hi | (Half(1) << (i - Half::BITS));
                    } else {
                        // Бит попадает в младшую половину
                        q.m_lo = q.m_lo | (Half(1) << i);
                    }
                }
            }
            current_den = current_den >> 1;
        }

        if constexpr (need_r) {
            if (rem_out)
                *rem_out = rem;
        }
        return q;
    }

    // --- Алгоритм D Кнута: деление "широкого" числа на "половинку" ---

    /**
     * @brief Длинное деление (алгоритм D Кнута) в системе счисления 2^LIMB_BITS.
     * @details Делит u (m лимбов) на v (n лимбов). v[n-1] != 0, n >= 1.
     *          u изменяется на месте. q получает m-n+1 лимбов (остальные
     *          должны быть предварительно обнулены), r — n лимбов.
     *          Разрядность лимба (64 или 32 бита) выбирается на этапе
     *          компиляции — так код работает и на MSVC без __int128.
     */
    static void divmnu(limb_t *q, limb_t *r, limb_t *u,
                       const limb_t *v, int m, int n) noexcept
    {
        // Делитель из одного лимба — простой проход.
        if (n == 1)
        {
            limb_t k = 0;
            for (int j = m - 1; j >= 0; --j)
            {
                wide_t cur = (static_cast<wide_t>(k) << LIMB_BITS) | u[j];
                q[j] = static_cast<limb_t>(cur / v[0]);
                k = static_cast<limb_t>(cur % v[0]);
            }
            r[0] = k;
            return;
        }

        constexpr wide_t B = static_cast<wide_t>(1) << LIMB_BITS;
        const int s = static_cast<int>(std::countl_zero(v[n - 1]));

        // Нормализация: сдвигаем делитель так, чтобы старший бит старшего лимба был установлен.
        limb_t vn[LIMBS];
        for (int i = n - 1; i > 0; --i)
            vn[i] = (v[i] << s) | (s ? (v[i - 1] >> (LIMB_BITS - s)) : 0);
        vn[0] = static_cast<limb_t>(v[0] << s);

        limb_t un[LIMBS + 1];
        un[m] = s ? (u[m - 1] >> (LIMB_BITS - s)) : 0;
        for (int i = m - 1; i > 0; --i)
            un[i] = (u[i] << s) | (s ? (u[i - 1] >> (LIMB_BITS - s)) : 0);
        un[0] = static_cast<limb_t>(u[0] << s);

        for (int j = m - n; j >= 0; --j)
        {
            // Оценка очередного лимба частного с коррекцией.
            wide_t num = (static_cast<wide_t>(un[j + n]) << LIMB_BITS) | un[j + n - 1];
            wide_t qhat = num / vn[n - 1];
            wide_t rhat = num % vn[n - 1];
            while (qhat >= B || qhat * vn[n - 2] > ((rhat << LIMB_BITS) | un[j + n - 2]))
            {
                --qhat;
                rhat += vn[n - 1];
                if (rhat >= B)
                    break;
            }

            // Умножение и вычитание: un[j .. j+n] -= qhat * vn
            swide_t k = 0;
            for (int i = 0; i < n; ++i)
            {
                wide_t p = qhat * vn[i];
                swide_t t = static_cast<swide_t>(un[i + j]) - k - static_cast<swide_t>(static_cast<limb_t>(p));
                un[i + j] = static_cast<limb_t>(t);
                k = static_cast<swide_t>(p >> LIMB_BITS) - (t >> LIMB_BITS);
            }
            swide_t t2 = static_cast<swide_t>(un[j + n]) - k;
            un[j + n] = static_cast<limb_t>(t2);
            q[j] = static_cast<limb_t>(qhat);

            if (t2 < 0)
            {
                // Оценка была на единицу больше — прибавляем делитель обратно.
                q[j] -= 1;
                limb_t carry = 0;
                for (int i = 0; i < n; ++i)
                {
                    wide_t ss = static_cast<wide_t>(un[i + j]) + vn[i] + carry;
                    un[i + j] = static_cast<limb_t>(ss);
                    carry = static_cast<limb_t>(ss >> LIMB_BITS);
                }
                un[j + n] = static_cast<limb_t>(un[j + n] + carry);
            }
        }

        // Денормализация остатка.
        for (int i = 0; i < n; ++i)
            r[i] = (un[i] >> s) | (s ? (un[i + 1] << (LIMB_BITS - s)) : 0);
    }

    struct DivResult;

    template <typename H>
    friend typename UBig<H>::DivResult operator/(const UBig<H>& a, const UBig<H>& b);

    template <typename H>
    friend UBig<H> operator%(const UBig<H>& a, const UBig<H>& b);

    UBig& operator/=(const UBig& o);
    UBig& operator%=(const UBig& o);

    static UBig fromString(const std::string &str)
    {
        UBig result(0);
        size_t idx = 0;

        // 1. Пропуск начальных пробелов
        while (idx < str.size()
               && (str[idx] == ' ' || str[idx] == '\t' || str[idx] == '\n' || str[idx] == '\r')) {
            idx++;
        }

        // 2. Пропуск знака '+'
        if (idx < str.size() && str[idx] == '+') {
            idx++;
        }

        // Степени десятки от 10^0 до 10^18 для аппаратного uint64_t
        static constexpr uint64_t POW10[] = {
            1ULL,
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
            1'000'000'000'000'000'000ULL // 10^18
        };

        while (idx < str.size()) {
            uint64_t chunk = 0;
            size_t count = 0;

            // Набираем до 18 цифр в обычный 64-битный регистр
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

            if (count == 0) {
                break; // Встретился нецифровой символ
            }

            // Одно длинное умножение и сложение на каждые 18 цифр
            result = result * UBig(POW10[count]) + UBig(chunk);
        }

        return result;
    }

    std::string toString() const
    {
        if (*this == UBig(0))
            return "0";

        // Быстрый выход для чисел, помещающихся в стандартный 64-битный int
        if (this->fits_u64()) {
            return std::to_string(this->to_u64());
        }

        static constexpr uint64_t BASE18 = 1'000'000'000'000'000'000ULL; // 10^18
        const UBig base_big(BASE18);

        std::string s;
        UBig cur = *this;

        while (cur > UBig(0)) {
            // Если осталась только младшая половина Half, передаем форматирование ей
            if (cur.high() == Half(0)) {
                std::string head = cur.low().toString();
                return head + s;
            }

            UBig rem(0);
            cur = divide<true, true>(cur, base_big, &rem);
            uint64_t chunk_val = rem.to_u64();

            // Переводим 64-битный остаток в строку
            std::string chunk = std::to_string(chunk_val);

            // Если частное еще не ноль, выравниваем блок ведущими нулями ровно до 18 знаков
            if (cur > UBig(0)) {
                if (chunk.size() < 18) {
                    chunk = std::string(18 - chunk.size(), '0') + chunk;
                }
            }

            s = chunk + s;
        }

        return s;
    }

    friend std::ostream &operator<<(std::ostream &os, const UBig &val)
    {
        return os << val.toString();
    }
};

template <typename Half>
struct UBig<Half>::DivResult {
    UBig<Half> first;
    UBig<Half> second;

    constexpr operator UBig<Half>() const { return first; }
};

template <typename Half>
inline typename UBig<Half>::DivResult operator/(const UBig<Half>& a, const UBig<Half>& b) {
    UBig<Half> r(0);
    UBig<Half> q = UBig<Half>::template divide<true, true>(a, b, &r);
    return typename UBig<Half>::DivResult{q, r};
}

template <typename Half>
inline UBig<Half> operator%(const UBig<Half>& a, const UBig<Half>& b) {
    UBig<Half> r(0);
    UBig<Half>::template divide<false, true>(a, b, &r);
    return r;
}

template <typename Half>
inline typename UBig<Half>::DivResult operator/(const UBig<Half>& a, const Half& b) {
    return a / UBig<Half>(b);
}

template <typename Half>
inline UBig<Half> operator%(const UBig<Half>& a, const Half& b) {
    return a % UBig<Half>(b);
}

template <typename Half>
inline constexpr bool operator<(const UBig<Half>& a, const Half& b) {
    return a.high() == Half(0) && a.low() < b;
}

template <typename Half>
inline constexpr bool operator<(const Half& a, const UBig<Half>& b) {
    return b.high() > Half(0) || a < b.low();
}

template <typename Half>
inline constexpr bool operator==(const UBig<Half>& a, const Half& b) {
    return a.high() == Half(0) && a.low() == b;
}

template <typename Half>
inline constexpr bool operator==(const Half& a, const UBig<Half>& b) {
    return b == a;
}

template <typename Half>
inline UBig<Half>& UBig<Half>::operator/=(const UBig<Half>& o) {
    return *this = (*this / o).first;
}

template <typename Half>
inline UBig<Half>& UBig<Half>::operator%=(const UBig<Half>& o) {
    return *this = *this % o;
}

using U256  = UBig<U128>;
using U512  = UBig<U256>;
using U1024 = UBig<U512>;
using U2048 = UBig<U1024>;
using U4096 = UBig<U2048>;

} // namespace bignum

using U256  = bignum::U256;
using U512  = bignum::U512;
using U1024 = bignum::U1024;
using U2048 = bignum::U2048;
using U4096 = bignum::U4096;
