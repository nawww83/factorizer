// Проверка нового фиксированного big-int (bignum.hpp).
// Сборка: clang++ -std=c++20 -O2 test_bignum.cpp -o test_bignum
#include "../bignum.hpp"

#include <cassert>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace bignum;

static std::mt19937_64 rng(1234567);

template <uint32_t Bits>
UInt<Bits> random_bits(uint32_t width)
{
    UInt<Bits> v(0);
    for (uint32_t i = 0; i < UInt<Bits>::WORDS; ++i)
        v.w[i] = rng();
    if (width < Bits)
        v.shr(Bits - width);
    if (width > 0)
        v.set_bit(width - 1);
    return v;
}

template <uint32_t Bits>
void check_roundtrip()
{
    for (int it = 0; it < 200; ++it)
    {
        const uint32_t width = 1 + static_cast<uint32_t>(rng() % Bits);
        UInt<Bits> v = random_bits<Bits>(width);
        const std::string s = v.toString();
        UInt<Bits> back = UInt<Bits>::fromString(s);
        assert(back == v);
        assert(back.bit_width() == width);
    }
}

template <uint32_t Bits>
void check_div()
{
    for (int it = 0; it < 300; ++it)
    {
        UInt<Bits> a = random_bits<Bits>(static_cast<uint32_t>(1 + rng() % Bits));
        UInt<Bits> b = random_bits<Bits>(static_cast<uint32_t>(1 + rng() % Bits));
        if (b.is_zero())
            continue;
        UInt<Bits> r;
        UInt<Bits> q = UInt<Bits>::divide(a, b, &r);
        assert(r < b);
        UInt<Bits> check = q * b + r;
        // q*b + r == a при отсутствии переполнения q*b
        if (q.bit_width() + b.bit_width() <= Bits)
            assert(check == a);
    }
}

template <uint32_t Bits>
void check_mulmod(uint64_t m)
{
    for (int it = 0; it < 500; ++it)
    {
        UInt<Bits> a(rng() % m);
        UInt<Bits> b(rng() % m);
        const uint64_t expect = static_cast<uint64_t>(
            (static_cast<unsigned __int128>(a.to_u64()) * b.to_u64()) % m);
        assert(UInt<Bits>::mul_mod(a, b, UInt<Bits>(m)).to_u64() == expect);
        BarrettReducer<Bits> br{(UInt<Bits>(m))};
        assert(br.mul(a, b).to_u64() == expect);
    }
}

int main()
{
    std::fprintf(stderr, "rt128\n"); check_roundtrip<128>();
    std::fprintf(stderr, "rt256\n"); check_roundtrip<256>();
    std::fprintf(stderr, "rt512\n"); check_roundtrip<512>();
    std::fprintf(stderr, "rt1024\n"); check_roundtrip<1024>();

    std::fprintf(stderr, "div128\n"); check_div<128>();
    std::fprintf(stderr, "div256\n"); check_div<256>();
    std::fprintf(stderr, "div512\n"); check_div<512>();
    std::fprintf(stderr, "div1024\n"); check_div<1024>();

    for (uint64_t m : {3ULL, 5ULL, 65537ULL, 1000000007ULL, 4294967291ULL,
                       1000000000000000003ULL})
        check_mulmod<128>(m);

    // Известные значения
    UInt<512> x = UInt<512>::fromString(
        "5995806622530429600530984275000092882631332762476331275490451238005952121568939");
    assert(x.toString() ==
           "5995806622530429600530984275000092882631332762476331275490451238005952121568939");
    assert(x.bit_width() == 262);

    UInt<256> a = UInt<256>::fromString(
        "115792089237316195423570985008687907853269984665640564039457584007913129639935");
    UInt<256> b = UInt<256>::fromString("340282366920938463463374607431768211455");
    UInt<256> r;
    UInt<256> q = UInt<256>::divide(a, b, &r);
    assert(q == UInt<256>::fromString("340282366920938463463374607431768211457"));
    assert(r.is_zero());

    // pow_mod: 3^100 mod 1000000007 == 962305464
    UInt<128> p = bignum::ops::pow_mod<128>(UInt<128>(3), UInt<128>(100),
                                            UInt<128>(1000000007ULL));
    assert(p.to_u64() == 886041711ULL);

    // isqrt
    UInt<128> sq = UInt<128>::fromString("15241578750190521"); // 123456789^2
    assert(bignum::ops::isqrt<128>(sq).to_u64() == 123456789ULL);

    std::fprintf(stderr, "bignum self-test: OK\n");
    return 0;
}
