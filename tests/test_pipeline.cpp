// Дымовой тест конвейера: factor_parallel_generic через решето siqs_xor.
// Проверяет, что новый SIQS встроен в общий конвейер (пробное деление,
// Поллард, решето) и что произведение найденных множителей равно исходному
// числу — и в авторежиме, и при принудительном выборе метода.
#include "generic_factorizer.hpp"

#include <chrono>
#include <iostream>
#include <string>

namespace {

int g_failures = 0;

template <typename T>
void run_case(const char *label, const std::string &str, int method, bool skip_ecm)
{
    T n = T::fromString(str);
    std::cout << "== " << label << " (" << n.bit_width() << " бит, method=" << method
              << ", skip_ecm=" << skip_ecm << ")\n";

    const auto t0 = std::chrono::steady_clock::now();
    std::string last_status;
    auto status = [&](const std::string &s) {
        if (s == last_status)
            return;
        last_status = s;
        std::cout << "   [status] " << s << "\n";
    };

    auto factors = bignum::generic_ecm::factor_parallel_generic(n, status, nullptr, skip_ecm, method);
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    T prod(1);
    for (const auto &f : factors) {
        std::cout << "   factor: " << f.toString() << "\n";
        prod = prod * f;
    }

    const bool ok = (prod == n) && factors.size() >= 2;
    std::cout << "   time: " << secs << " s, decomposition ok: " << ok << "\n\n";
    if (!ok)
        ++g_failures;
}

} // namespace

int main()
{
    // Полупростое 562949953433701 * 1125899906897009.
    const std::string semiprime100 = "633825300158679507911716700309";
    // 205-битное полупростое из эталонного прогона tests/siqs_xor.cpp.
    const std::string semiprime205 =
        "33801295533024112192430315909773114733482083162238224697211891";

    run_case<U128>("100-бит, авто, без ECM", semiprime100, bignum::generic_ecm::kSieveAuto, true);
    run_case<U128>("100-бит, принудительный SIQS", semiprime100,
                   bignum::generic_ecm::kSieveForced, false);
    run_case<U256>("205-бит, авто, без ECM", semiprime205, bignum::generic_ecm::kSieveAuto, true);

    std::cout << (g_failures == 0 ? "pipeline self-test: OK\n" : "pipeline self-test: FAILED\n");
    return g_failures == 0 ? 0 : 1;
}
