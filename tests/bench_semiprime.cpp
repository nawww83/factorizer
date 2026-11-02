// Бенчмарк решета.
//
//   bench_semiprime <бит> [<бит> ...]     — сгенерировать полупростые p×q
//                                           заданной битности и разложить их;
//   bench_semiprime --n <десятичное>      — разложить заданное число.
//
// Факторизация идёт тем же конвейером, что и в приложении: режим «Авто»,
// ECM пропущен (типичный сценарий для полупростых чисел).
#include "generic_factorizer.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#if defined(_WIN32)
#  include <windows.h>
#  include <psapi.h>
#endif

namespace {

// Пиковый рабочий набор текущего процесса (МБ); 0, если ОС не дала данных.
double peak_working_set_mb()
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return static_cast<double>(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0);
#endif
    return 0.0;
}

template <typename T>
int measure(const T &n)
{
    const auto t0 = std::chrono::steady_clock::now();
    auto factors = bignum::generic_ecm::factor_parallel_generic(
        n, nullptr, nullptr, /*skip_ecm=*/true, bignum::generic_ecm::kSieveAuto);
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    T prod(1);
    for (const auto &f : factors)
        prod = prod * f;

    std::cout << "bits=" << n.bit_width()
              << " time=" << secs
              << " peakMB=" << peak_working_set_mb()
              << " factors=" << factors.size()
              << " ok=" << (prod == n) << "\n  n=" << n.toString() << "\n";
    for (const auto &f : factors)
        std::cout << "  factor=" << f.toString() << "\n";
    return prod == n ? 0 : 3;
}

template <typename T>
int run_semiprime(int bits, uint64_t seed)
{
    T p(0), q(0);
    const T n = bignum::generic_ecm::generate_semiprime<T>(static_cast<unsigned>(bits),
                                                           &p, &q, nullptr, seed);
    if (n == T(0)) {
        std::cerr << "не удалось сгенерировать полупростое " << bits << " бит\n";
        return 2;
    }
    std::cout << "  p=" << p.toString() << "\n  q=" << q.toString() << "\n";
    return measure(n);
}

template <typename T>
int run_number(const std::string &dec)
{
    return measure(T::fromString(dec));
}

int run_any_number(const std::string &dec)
{
    // Тип выбираем по фактической битности числа (предел решета — 512 бит).
    const U4096 probe = U4096::fromString(dec);
    const uint32_t bits = probe.bit_width();
    if (bits <= 128)
        return run_number<U128>(dec);
    if (bits <= 256)
        return run_number<U256>(dec);
    if (bits <= 512)
        return run_number<U512>(dec);
    std::cerr << "число длиннее 512 бит: решето неприменимо\n";
    return 1;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "usage: bench_semiprime <bits> [<bits> ...] | --n <число>\n";
        return 1;
    }

    if (std::string(argv[1]) == "--n") {
        if (argc < 3) {
            std::cerr << "usage: bench_semiprime --n <число>\n";
            return 1;
        }
        return run_any_number(argv[2]);
    }

    for (int i = 1; i < argc; ++i) {
        const int bits = std::atoi(argv[i]);
        if (bits < 8 || bits > 512) {
            std::cerr << "битность должна быть в диапазоне 8..512\n";
            return 1;
        }
        const uint64_t seed = 1000ULL + static_cast<uint64_t>(bits);

        std::cout << "== " << bits << " бит\n";
        int rc = 1;
        if (bits <= 128)
            rc = run_semiprime<U128>(bits, seed);
        else if (bits <= 256)
            rc = run_semiprime<U256>(bits, seed);
        else
            rc = run_semiprime<U512>(bits, seed);
        if (rc != 0)
            return rc;
    }
    return 0;
}
