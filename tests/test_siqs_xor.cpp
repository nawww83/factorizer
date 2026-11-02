// Проверка нового решета: сборка clang++ -std=c++20 -O2 test_siqs_xor.cpp
#include "../siqs_xor.hpp"

#include <cassert>
#include <cstdio>
#include <string>

using namespace bignum;

static bool check(const std::string &dec)
{
    UInt<512> n = UInt<512>::fromString(dec);
    std::printf("Факторизация %u бит...\n", n.bit_width());
    auto d = siqs_xor::factor<UInt<512>>(n, [](const std::string &s)
                                         { std::printf("  %s\n", s.c_str()); });
    if (!d.has_value())
    {
        std::printf("  делитель не найден\n");
        return false;
    }
    UInt<512> q;
    UInt<512>::divide(n, *d, &q);
    std::printf("  p = %s\n  q = %s\n", d->toString().c_str(), q.toString().c_str());
    return *d > UInt<512>(1) && *d < n && (n % *d).is_zero();
}

int main()
{
    ::u128::Globals::SetStop(false);

    // 205-битное полупростое из эталонного прогона.
    const bool ok1 = check("33801295533024112192430315909773114733482083162238224697211891");
    assert(ok1);

    // 100-битное полупростое (два 50-битных простых).
    const bool ok2 = check("633825300158679507911716700309");
    assert(ok2);

    std::puts("siqs_xor self-test: OK");
    return 0;
}
