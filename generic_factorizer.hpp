#pragma once

#include "bignum.hpp"
#include "siqs_xor.hpp"
#include "u128_utils.h"
#include "cpu_topology.hpp"

#include <vector>
#include <string>
#include <functional>
#include <atomic>
#include <mutex>
#include <thread>
#include <optional>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <utility>
#include <random>

namespace bignum {

struct StageProfiler {
    using Clock = std::chrono::steady_clock;
    Clock::time_point start_time = Clock::now();
    Clock::time_point stage_start = Clock::now();
    std::string current_stage;
    std::vector<std::pair<std::string, double>> stage_durations;

    void start_stage(const std::string& name) {
        auto now = Clock::now();
        if (!current_stage.empty()) {
            double elapsed = std::chrono::duration<double>(now - stage_start).count();
            record_elapsed(current_stage, elapsed);
        }
        current_stage = name;
        stage_start = now;
    }

    void finish() {
        auto now = Clock::now();
        if (!current_stage.empty()) {
            double elapsed = std::chrono::duration<double>(now - stage_start).count();
            record_elapsed(current_stage, elapsed);
            current_stage.clear();
        }
    }

    double total_time() const {
        return std::chrono::duration<double>(Clock::now() - start_time).count();
    }

    std::string format_report() const {
        if (stage_durations.empty()) return "";
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(3);
        ss << "\n--- Затраченное время по шагам ---\n";
        for (const auto& [stage, duration] : stage_durations) {
            ss << "  • " << stage << ": " << duration << " с\n";
        }
        return ss.str();
    }

private:
    void record_elapsed(const std::string& name, double elapsed) {
        for (auto& [stage, duration] : stage_durations) {
            if (stage == name) {
                duration += elapsed;
                return;
            }
        }
        stage_durations.emplace_back(name, elapsed);
    }
};

} // namespace bignum

namespace bignum::generic_ecm {

// --- Вспомогательные операции сложения/вычитания по модулю N ---

template <typename T>
inline T mod_add(const T& a, const T& b, const T& m) {
    T res = a + b;
    if (res >= m || res < a) {
        res = res - m;
    }
    return res;
}

template <typename T>
inline T mod_sub(const T& a, const T& b, const T& m) {
    if (a >= b) {
        return a - b;
    } else {
        return m - (b - a);
    }
}

template <typename T>
bool miller_test_generic(T d, const T& n) {
    static thread_local std::mt19937_64 rng(1337);
    T a = T(rng()) ^ (T(rng()) << 64);
    a = ::bignum::ops::div_r(a, n - T(4)) + T(2); // a in [2, n - 2]

    T x = a;
    ::bignum::ops::int_power_mod(x, d, n);

    if (x == T(1) || x == (n - T(1))) {
        return true;
    }

    T n_minus_1 = n - T(1);
    while (d != n_minus_1) {
        ::bignum::ops::square_mod(x, n);
        d = d << 1;

        if (x == T(1)) {
            return false;
        }
        if (x == n_minus_1) {
            return true;
        }
    }
    return false;
}

template<typename T>
bool is_prime_generic(const T &n, int k = 20)
{
    if (n <= T(1)) return false;
    if (n <= T(3)) return true;
    if ((n & T(1)) == T(0)) return false;

    T d = n - T(1);
    while ((d & T(1)) == T(0)) {
        d = d >> 1;
    }

    for (int i = 0; i < k; ++i) {
        if (!miller_test_generic(d, n)) {
            return false;
        }
    }
    return true;
}

// --- 0. Pollard-Brent rho с пакетированием GCD ---

template <typename T>
T pollard_brent_rho(const T& n, uint64_t max_steps = 200000) {
    if (n <= T(1)) return n;
    if ((n & T(1)) == T(0)) return T(2);

    T y = T(2);
    T c = T(1);
    T g = T(1);
    T q = T(1);
    T ys = T(0);
    T x = T(0);

    uint64_t r = 1;
    uint64_t steps = 0;

    while (g == T(1) && steps < max_steps && !::u128::Globals::LoadStop()) {
        x = y;
        // r растёт как 2^k и на последних витках достигает ~100000, поэтому этот
        // цикл — самый долгий неотменяемый участок ранней стадии. Проверяем
        // стоп-флаг прямо в нём, иначе кнопка «Остановить» на крупных числах
        // «залипает» на секунды.
        for (uint64_t i = 0; i < r; ++i) {
            if ((i & 0x3FF) == 0 && ::u128::Globals::LoadStop())
                return T(1);
            ::bignum::ops::square_mod(y, n);
            y = mod_add(y, c, n);
        }

        uint64_t k = 0;
        while (k < r && g == T(1) && steps < max_steps && !::u128::Globals::LoadStop()) {
            ys = y;
            uint64_t limit = std::min(static_cast<uint64_t>(128), r - k);
            for (uint64_t i = 0; i < limit; ++i) {
                ::bignum::ops::square_mod(y, n);
                y = mod_add(y, c, n);
                T diff = (x >= y) ? (x - y) : (y - x);
                ::bignum::ops::mult_mod(q, diff, n);
            }

            g = ::bignum::ops::gcd(q, n);
            k += limit;
            steps += limit;
        }

        r *= 2;
    }

    if (g == n) {
        do {
            if (::u128::Globals::LoadStop())
                return T(1);
            ::bignum::ops::square_mod(ys, n);
            ys = mod_add(ys, c, n);
            T diff = (x >= ys) ? (x - ys) : (ys - x);
            g = ::bignum::ops::gcd(diff, n);
        } while (g == T(1));
    }

    if (g > T(1) && g < n) {
        return g;
    }

    return T(1);
}

// --- 1. Pollard p - 1 ---

template <typename T>
std::optional<T> pollard_pm1(const T& n, uint32_t B1 = 100000) {
    auto primes = ::u128::utils::primes(B1);
    T a = T(2);

    for (uint32_t p : primes) {
        if (::u128::Globals::LoadStop()) return std::nullopt;

        uint32_t p_pow = p;
        while (1ULL * p_pow * p <= B1) {
            p_pow *= p;
        }

        ::bignum::ops::int_power_mod(a, T(p_pow), n);

        T diff = (a > T(1)) ? (a - T(1)) : (n - T(1));
        T g = ::bignum::ops::gcd(diff, n);
        if (g > T(1) && g < n) {
            return g;
        }
        if (g == n) {
            break;
        }
    }
    return std::nullopt;
}

// --- 2. ECM: Арифметика кривых Монтгомери ---

template <typename T>
struct MontgomeryPoint {
    T X;
    T Z;
};

// --- Безопасное удвоение по модулю N без переполнения ---
template <typename T>
inline T mod_double(const T& a, const T& m) {
    return mod_add(a, a, m);
}

// --- Безопасный бинарный инвертор a^(-1) mod m (для нечетного m) ---
template <typename T>
std::optional<T> inv_mod_binary(const T& a_in, const T& m) {
    if (m <= T(1)) return std::nullopt;
    T a = ::bignum::ops::div_r(a_in, m);
    if (a == T(0)) return std::nullopt;

    T u = a;
    T v = m;
    T x1 = T(1);
    T x2 = T(0);

    // Вспомогательная лямбда: безопасно вычисляет (x + m) >> 1 без переполнения разрядности
    auto safe_add_m_div2 = [](const T& x, const T& mod) -> T {
        // Так как x нечетное и mod нечетное: (x + mod) / 2 = (x >> 1) + (mod >> 1) + 1
        return (x >> 1) + (mod >> 1) + T(1);
    };

    while (u != T(0) && v != T(0)) {
        while ((u & T(1)) == T(0)) {
            u = u >> 1;
            if ((x1 & T(1)) != T(0)) {
                x1 = safe_add_m_div2(x1, m);
            } else {
                x1 = x1 >> 1;
            }
        }

        while ((v & T(1)) == T(0)) {
            v = v >> 1;
            if ((x2 & T(1)) != T(0)) {
                x2 = safe_add_m_div2(x2, m);
            } else {
                x2 = x2 >> 1;
            }
        }

        if (u >= v) {
            u = u - v;
            x1 = mod_sub(x1, x2, m);
        } else {
            v = v - u;
            x2 = mod_sub(x2, x1, m);
        }
    }

    if (v == T(1)) return x2;
    if (u == T(1)) return x1;

    // НОД(a, m) > 1: делитель перехватывается снаружи через gcd
    return std::nullopt;
}

// --- Канонические формулы Монтгомери (при условии Z_diff = 1) ---

template <typename T>
inline MontgomeryPoint<T> montgomery_double(const MontgomeryPoint<T>& P,
                                            const T& A24,
                                            const BarrettReducer<T::BITS>& br) {
    const T& n = br.modulus();
    T u = mod_add(P.X, P.Z, n);
    T v = mod_sub(P.X, P.Z, n);

    T u2 = br.sqr(u);
    T v2 = br.sqr(v);
    T diff = mod_sub(u2, v2, n); // 4 * X * Z

    T X2 = br.mul(u2, v2);

    T term = br.mul(diff, A24);
    term = mod_add(v2, term, n);

    T Z2 = br.mul(diff, term);

    return {X2, Z2};
}

template <typename T>
inline MontgomeryPoint<T> montgomery_add(const MontgomeryPoint<T>& P,
                                         const MontgomeryPoint<T>& Q,
                                         const T& x_diff,
                                         const BarrettReducer<T::BITS>& br) {
    const T& n = br.modulus();
    T u1 = mod_add(P.X, P.Z, n);
    T v1 = mod_sub(P.X, P.Z, n);
    T u2 = mod_add(Q.X, Q.Z, n);
    T v2 = mod_sub(Q.X, Q.Z, n);

    T t1 = br.mul(u1, v2);
    T t2 = br.mul(v1, u2);

    T add_t = mod_add(t1, t2, n);
    T sub_t = mod_sub(t1, t2, n);

    T X_out = br.sqr(add_t);
    T Z_out = br.sqr(sub_t);
    Z_out = br.mul(Z_out, x_diff);

    return {X_out, Z_out};
}


// --- Решето простых для Stage 2 ---
inline std::vector<uint8_t> make_prime_sieve(uint32_t max_val) {
    std::vector<uint8_t> is_prime(static_cast<size_t>(max_val) + 1, 1);
    is_prime[0] = 0;
    if (max_val >= 1) is_prime[1] = 0;
    for (uint32_t p = 2; 1ULL * p * p <= max_val; ++p) {
        if (is_prime[p]) {
            // 64-битный индекс исключает переполнение i += p у верхней границы диапазона.
            for (uint64_t i = 1ULL * p * p; i <= max_val; i += p) {
                is_prime[static_cast<size_t>(i)] = 0;
            }
        }
    }
    return is_prime;
}

// --- Умножение точки на скаляр через лестницу Монтгомери ---
template <typename T>
MontgomeryPoint<T> montgomery_multiply(uint32_t k,
                                       const MontgomeryPoint<T>& P,
                                       const T& A24,
                                       const BarrettReducer<T::BITS>& br) {
    MontgomeryPoint<T> R0 = P;
    MontgomeryPoint<T> R1 = montgomery_double(P, A24, br);

    int msb = 31 - std::countl_zero(k);
    for (int bit = msb - 1; bit >= 0; --bit) {
        if ((k >> bit) & 1) {
            R0 = montgomery_add(R1, R0, P.X, br);
            R1 = montgomery_double(R1, A24, br);
        } else {
            R1 = montgomery_add(R0, R1, P.X, br);
            R0 = montgomery_double(R0, A24, br);
        }
    }
    return R0;
}

// --- ECM Stage 2: Baby-Step Giant-Step (BSGS) с редуктором Барретта ---
template <typename T>
std::optional<T> ecm_stage2(MontgomeryPoint<T>& Q,
                            const T& A24,
                            const BarrettReducer<T::BITS>& br,
                            uint32_t B1,
                            uint32_t B2,
                            const std::vector<uint8_t>& is_prime_table) {
    const T& n = br.modulus();

    static const uint32_t J[] = {
        1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
        53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103
    };

    // 1. Baby-steps: вычисляем S_j = j * Q (в проективных координатах)
    MontgomeryPoint<T> baby[106];
    baby[1] = Q;

    for (uint32_t j : J) {
        if (j == 1) continue;
        baby[j] = montgomery_multiply(j, Q, A24, br);
        if (baby[j].Z == T(0)) {
            T g = ::bignum::ops::gcd(baby[j].X, n);
            if (g > T(1) && g < n) return g;
            return std::nullopt;
        }
    }

    // 2. Giant-steps setup: G1 = 210 * Q, G2 = 2 * G1
    MontgomeryPoint<T> G1 = montgomery_multiply(210, Q, A24, br);
    if (G1.Z == T(0)) return std::nullopt;

    auto g1_z_inv = inv_mod_binary(G1.Z, n);
    if (!g1_z_inv.has_value()) {
        T g = ::bignum::ops::gcd(G1.Z, n);
        if (g > T(1) && g < n) return g;
        return std::nullopt;
    }
    G1.X = br.mul(G1.X, *g1_z_inv);
    G1.Z = T(1);

    MontgomeryPoint<T> G_prev = G1;
    MontgomeryPoint<T> G_curr = montgomery_double(G1, A24, br);

    uint32_t m_max = (B2 + 105) / 210 + 1;
    T acc = T(1);
    int acc_count = 0;

    for (uint32_t m = 2; m <= m_max; ++m) {
        if (::u128::Globals::LoadStop()) return std::nullopt;

        if (m * 210 + 105 > B1) {
            for (uint32_t j : J) {
                bool need = false;
                uint32_t q_plus = m * 210 + j;
                if (q_plus > B1 && q_plus <= B2 && is_prime_table[q_plus]) need = true;
                if (!need && m * 210 >= j) {
                    uint32_t q_minus = m * 210 - j;
                    if (q_minus > B1 && q_minus <= B2 && is_prime_table[q_minus]) need = true;
                }

                if (need) {
                    // Перекрестное умножение без делений через Барретта:
                    T t1 = br.mul(G_curr.X, baby[j].Z);
                    T t2 = br.mul(baby[j].X, G_curr.Z);
                    T diff = (t1 >= t2) ? (t1 - t2) : (n - (t2 - t1));

                    if (diff != T(0)) {
                        acc = br.mul(acc, diff);
                        acc_count++;
                    }
                }
            }
        }

        if (acc_count >= 128) {
            T g = ::bignum::ops::gcd(acc, n);
            if (g > T(1) && g < n) return g;
            if (g == n) return std::nullopt;
            acc = T(1);
            acc_count = 0;
        }

        if (m < m_max) {
            MontgomeryPoint<T> G_next = montgomery_add(G_curr, G1, G_prev.X, br);
            G_prev = G_curr;
            G_curr = G_next;
        }
    }

    if (acc_count > 0) {
        T g = ::bignum::ops::gcd(acc, n);
        if (g > T(1) && g < n) return g;
    }

    return std::nullopt;
}

// --- Запуск одной кривой ECM (Stage 1 + Stage 2) ---
template <typename T>
std::optional<T> ecm_one_curve(const T& n, uint64_t seed, uint32_t B1, uint32_t B2,
                               const std::vector<unsigned>& primes,
                               const std::vector<uint8_t>& is_prime_table) {
    BarrettReducer<T::BITS> br(n); // Редуктор Барретта создается один раз на кривую

    uint64_t sigma_val = (seed % 1000000000ULL) + 6;
    T sigma(sigma_val);

    T sigma2 = br.sqr(sigma);
    T u = mod_sub(sigma2, T(5), n);
    T v = mod_double(mod_double(sigma, n), n);

    T u3 = br.mul(br.sqr(u), u);
    T v3 = br.mul(br.sqr(v), v);

    auto v3_inv = inv_mod_binary(v3, n);
    if (!v3_inv.has_value()) {
        T g = ::bignum::ops::gcd(v3, n);
        if (g > T(1) && g < n) return g;
        return std::nullopt;
    }
    T x0 = br.mul(u3, *v3_inv);
    MontgomeryPoint<T> P{x0, T(1)};

    T diff_vu = mod_sub(v, u, n);
    T diff_vu3 = br.mul(br.sqr(diff_vu), diff_vu);

    T sum_3uv = mod_add(mod_add(mod_double(u, n), u, n), v, n);
    T A_num = br.mul(diff_vu3, sum_3uv);

    T denom = br.mul(u3, v);
    denom = mod_double(mod_double(mod_double(mod_double(denom, n), n), n), n); // * 16

    auto denom_inv = inv_mod_binary(denom, n);
    if (!denom_inv.has_value()) {
        T g = ::bignum::ops::gcd(denom, n);
        if (g > T(1) && g < n) return g;
        return std::nullopt;
    }
    T A24 = br.mul(A_num, *denom_inv);

    // === Stage 1 с пакетированием скаляров ===
    uint64_t current_k = 1;

    for (unsigned p : primes) {
        if (p > B1) break;
        if (::u128::Globals::LoadStop()) return std::nullopt;

        uint32_t p_pow = p;
        while (1ULL * p_pow * p <= B1) {
            p_pow *= p;
        }

        if (current_k * p_pow > 0xFFFFFFFFULL) {
            P = montgomery_multiply(static_cast<uint32_t>(current_k), P, A24, br);
            if (P.Z == T(0)) return std::nullopt;

            auto z_inv = inv_mod_binary(P.Z, n);
            if (!z_inv.has_value()) {
                T g = ::bignum::ops::gcd(P.Z, n);
                if (g > T(1) && g < n) return g;
                return std::nullopt;
            }
            P.X = br.mul(P.X, *z_inv);
            P.Z = T(1);

            current_k = p_pow;
        } else {
            current_k *= p_pow;
        }
    }

    if (current_k > 1) {
        P = montgomery_multiply(static_cast<uint32_t>(current_k), P, A24, br);
        if (P.Z == T(0)) return std::nullopt;

        auto z_inv = inv_mod_binary(P.Z, n);
        if (!z_inv.has_value()) {
            T g = ::bignum::ops::gcd(P.Z, n);
            if (g > T(1) && g < n) return g;
            return std::nullopt;
        }
        P.X = br.mul(P.X, *z_inv);
        P.Z = T(1);
    }

    // === Stage 2 (BSGS) ===
    if (B2 > B1) {
        return ecm_stage2(P, A24, br, B1, B2, is_prime_table);
    }

    return std::nullopt;
}


// --- Параллельный диспетчер ECM ---
template <typename T>
std::optional<T> ecm_parallel(const T& n,
                              unsigned max_curves,
                              uint32_t B1,
                              uint32_t B2,
                              std::function<void(const std::string&)> status_cb) {
    auto primes = ::u128::utils::primes(B1);
    auto is_prime_table = make_prime_sieve(B2);

    // Только P-ядра: E-ядра на гибридных процессорах заметно медленнее.
    // Потоки дополнительно закрепляются за P-ядрами (см. pin_current_thread).
    unsigned num_threads = ::bignum::cpu::worker_thread_count();

    std::atomic<bool> found{false};
    std::atomic<uint32_t> curve_counter{0};
    std::mutex res_mutex;
    std::optional<T> result = std::nullopt;

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (unsigned t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            ::bignum::cpu::pin_current_thread(t);
            while (!found.load(std::memory_order_relaxed) &&
                   !::u128::Globals::LoadStop()) {
                uint32_t c_idx = curve_counter.fetch_add(1);
                if (c_idx >= max_curves) break;

                if (status_cb && c_idx % 4 == 0) {
                    status_cb("ECM: Тестирование кривых (" + std::to_string(c_idx) +
                              " / " + std::to_string(max_curves) + "), B1=" + std::to_string(B1) +
                              ", B2=" + std::to_string(B2));
                }

                uint64_t seed = 0x9E3779B97F4A7C15ULL * (c_idx + 1) + t * 9973;
                auto factor = ecm_one_curve(n, seed, B1, B2, primes, is_prime_table);

                if (factor.has_value()) {
                    std::lock_guard<std::mutex> lock(res_mutex);
                    if (!result.has_value()) {
                        result = factor;
                        found.store(true, std::memory_order_relaxed);
                    }
                    return;
                }
            }
        });
    }

    for (auto& w : workers) {
        if (w.joinable()) w.join();
    }

    return result;
}

// --- 4. Главный конвейер факторизации ---
// ECM-лестница запускается только там, где решето неприменимо (в авторежиме —
// для чисел крупнее kSieveMaxBits). Перед решетом ECM не работает: на 205-битном
// полупростом лестница тратила ~85 с и делитель не находила, тогда как решето
// укладывалось в 12 с.
//
// Решето — это новый 2LP-SIQS (siqs_xor.hpp). Замеры на полупростых p×q:
//  - 100 бит: ~3 с; 205 бит: ~4 с на 8 потоках; 262 бита: единицы минут.
//  - Реализация работает на 512-битной арифметике, поэтому числа длиннее
//    512 бит ею не раскладываются (для них остаются ECM и методы Полларда).
//  - В авторежиме решето запускается до kSieveMaxBits (320 бит): выше время
//    растёт очень быстро, и сутки на число — не результат.
constexpr unsigned kSieveMaxBits = 320;      // авторежим
constexpr unsigned kSieveHardMaxBits = 512;  // жёсткий предел реализации
constexpr unsigned kSiqsMaxBits = kSieveMaxBits; // совместимость с UI-сообщениями

// Способ выбора решета (значение совпадает с индексом в GUI):
//   0 — авто: решето до kSieveMaxBits, ECM-лестница для чисел крупнее;
//   1 — принудительно решето: без ограничения по разрядности (пользователь сам
//       решает, что готов ждать; выше 512 бит метод неприменим).
enum SieveChoice : int {
    kSieveAuto = 0,
    kSieveForced = 1,
    // Устаревшее имя, оставлено для совместимости исходников.
    kSieveSiqs = kSieveForced,
};

template <typename T>
std::vector<T> factor_parallel_generic(T n,
                                       std::function<void(const std::string&)> status_cb = nullptr,
                                       ::bignum::StageProfiler* profiler = nullptr,
                                       int sieve_method = kSieveAuto,
                                       std::function<void(double, const std::string&)> progress_cb = nullptr) {
    std::vector<T> factors;

    // Единая точка отчёта о прогрессе: доля 0..1 от всей факторизации + фаза.
    // Вызывается из рабочих потоков; получатель обязан быть потокобезопасным.
    auto progress = [&](double fraction, const std::string& stage) {
        if (progress_cb)
            progress_cb(fraction, stage);
    };

    if (n <= T(1)) return factors;

    progress(0.01, "подготовка");

    // 0. Мгновенное снятие степеней двойки битовыми сдвигами
    if ((n & T(1)) == T(0)) {
        if (profiler) profiler->start_stage("Снятие четности");
        while ((n & T(1)) == T(0)) {
            factors.push_back(T(2));
            n = n >> 1;
        }
    }

    if (n == T(1)) {
        if (profiler) profiler->finish();
        return factors;
    }

    // generic_factorizer.hpp -> Шаг 1 (Пробное деление)
    if (status_cb)
        status_cb("Пробное деление...");
    if (profiler)
        profiler->start_stage("Пробное деление (до 65536)");
    progress(0.02, "пробное деление");

    auto small_primes = ::u128::utils::primes(65536);
    size_t prime_idx = 0;
    for (unsigned p : small_primes) {
        if (p == 2)
            continue;
        if (::u128::Globals::LoadStop())
            return factors;
        // Прогресс пробного деления: 2%..6% по доле перебранных простых.
        if ((++prime_idx & 1023) == 0 && !small_primes.empty())
            progress(0.02 + 0.04 * (static_cast<double>(prime_idx)
                                    / static_cast<double>(small_primes.size())),
                     "пробное деление");

        // Безопасное сравнение без умножения в типе T:
        // если n уложилось в uint64_t и p * p > n, дальше простых делителей быть не может
        if (n <= T(std::numeric_limits<uint64_t>::max())) {
            uint64_t n64 = ::u128::utils::to_u64_low(n);
            if (1ULL * p * p > n64)
                break;
        }

        if (::bignum::ops::rem_u32(n, p) == 0) {
            T P_val(p);
            do {
                factors.push_back(P_val);
                n = ::bignum::ops::div_by_u32(n, p);
            } while (::bignum::ops::rem_u32(n, p) == 0);

            if (n == T(1))
                break;
        }
    }

    if (n == T(1)) {
        if (profiler) profiler->finish();
        return factors;
    }

    // 2. Проверка на простоту перед тяжелыми методами
    progress(0.06, "проверка простоты");
    if (is_prime_generic(n)) {
        factors.push_back(n);
        if (profiler) profiler->finish();
        return factors;
    }

    // 3. Pollard p - 1 (быстрый этап)
    if (status_cb) status_cb("Pollard p - 1...");
    if (profiler) profiler->start_stage("Pollard p - 1");
    progress(0.08, "Pollard p - 1");
    auto pm1_res = pollard_pm1(n, 50000);
    if (pm1_res.has_value()) {
        T f = *pm1_res;
        auto sub1 = factor_parallel_generic(f, status_cb, profiler, sieve_method, progress_cb);
        auto sub2 = factor_parallel_generic(::bignum::ops::div_q(n, f), status_cb, profiler, sieve_method, progress_cb);
        factors.insert(factors.end(), sub1.begin(), sub1.end());
        factors.insert(factors.end(), sub2.begin(), sub2.end());
        return factors;
    }

    // 4. Pollard-Brent rho (для делителей < 45 бит)
    if (status_cb) status_cb("Pollard-Brent rho...");
    if (profiler) profiler->start_stage("Pollard-Brent rho");
    progress(0.11, "Pollard-Brent rho");
    T f_rho = ::bignum::generic_ecm::pollard_brent_rho(n, 200000);
    if (f_rho > T(1) && f_rho < n) {
        auto sub1 = factor_parallel_generic(f_rho, status_cb, profiler, sieve_method, progress_cb);
        auto sub2 = factor_parallel_generic(::bignum::ops::div_q(n, f_rho), status_cb, profiler, sieve_method, progress_cb);
        factors.insert(factors.end(), sub1.begin(), sub1.end());
        factors.insert(factors.end(), sub2.begin(), sub2.end());
        return factors;
    }

    // 5. Решето как единая точка входа: пробуем 2LP-SIQS (siqs_xor.hpp).
    //    Возвращает true, если делитель найден и обе части уже разложены —
    //    тогда факторизация завершена. Проверка предела разрядности здесь же:
    //    реализация работает на 512-битной арифметике.
    const uint32_t n_bits = n.bit_width();
    auto try_sieve = [&](bool force) -> bool {
        if (n_bits > kSieveHardMaxBits) {
            if (status_cb)
                status_cb("SIQS пропущен: " + std::to_string(n_bits) + " бит превышает предел "
                          + std::to_string(kSieveHardMaxBits) + " бит.");
            return false;
        }

        if (status_cb)
            status_cb(force ? "Запуск SIQS (принудительно, без ограничения разрядности)..."
                            : "Запуск SIQS (квадратичное решето)...");
        if (profiler)
            profiler->start_stage("SIQS (Квадратичное решето)");

        // Локальный прогресс решета (0..1) отображается в полосу 14%..96%.
        auto sieve_progress = [&](double local, const std::string &stage) {
            progress(0.14 + 0.82 * local, stage);
        };
        auto sieve_res = ::bignum::siqs_xor::factor(n, status_cb, sieve_progress);
        if (!sieve_res.has_value())
            return false;

        T f = *sieve_res;
        if (f <= T(1) || f >= n)
            return false;

        auto sub1 = factor_parallel_generic(f, status_cb, profiler, sieve_method, progress_cb);
        auto sub2 = factor_parallel_generic(::bignum::ops::div_q(n, f), status_cb, profiler, sieve_method, progress_cb);
        factors.insert(factors.end(), sub1.begin(), sub1.end());
        factors.insert(factors.end(), sub2.begin(), sub2.end());
        return true;
    };

    // Быстрый прогон отдельным шагом не нужен: авторежим и так доходит до решета
    // ниже (п. 8), а повторный вызов на тех же данных с фиксированным сидом давал
    // бы тот же (детерминированный) результат.

    // 6. Полноценная лестница двухэтапного ECM (Stage 1 + Stage 2) — только для
    // чисел крупнее kSieveMaxBits, где решето в авторежиме всё равно не запустится.
    // Перед решетом ECM не работает: для полупростых p×q она не находит делитель,
    // и её время целиком уходит в накладные расходы (см. замеры в README).
    // При явном выборе решета лестница ECM тоже не запускается: пользователь
    // просил конкретный метод, а на 512 битах ECM-план верхних тиров идёт часами
    // до собственно решета.
    if (sieve_method == kSieveAuto && n.bit_width() > kSieveMaxBits) {
        struct EcmStage
        {
            unsigned curves;
            uint32_t B1;
            uint32_t B2;
        };

        size_t bits = n.bit_width();
        std::vector<EcmStage> ecm_plan;

        // Обязательная часть плана (битность здесь всегда > 320): сначала дешёвые
        // тиры для мелких факторов, затем средние — 150 и 450 кривых.
        ecm_plan.push_back({25, 2000, 60000});
        ecm_plan.push_back({50, 10000, 300000});
        ecm_plan.push_back({150, 50000, 1500000});
        ecm_plan.push_back({450, 110000, 3500000});

        // Продвинутый план для больших чисел (> 350 бит)
        if (bits > 350) {
            ecm_plan.push_back({1000, 250000, 8500000});   // Нацелен на поиск факторов ~35-40 бит
            ecm_plan.push_back({3000, 1000000, 40000000}); // Нацелен на поиск факторов ~45-50 бит
        }

        // Тяжелый план для экстремальных чисел (> 600 бит)
        if (bits > 600) {
            ecm_plan.push_back({10000, 3000000, 120000000});  // Факторы ~55-60 бит
            ecm_plan.push_back({25000, 11000000, 500000000}); // Факторы ~65-70 бит
        }

        for (const auto &stage : ecm_plan) {
            if (::u128::Globals::LoadStop())
                return factors;
            std::string stage_name = "ECM (B1=" + std::to_string(stage.B1)
                                     + ", B2=" + std::to_string(stage.B2)
                                     + ", кривых=" + std::to_string(stage.curves) + ")";
            if (profiler)
                profiler->start_stage(stage_name);

            progress(0.14 + 0.82 * (static_cast<double>(&stage - ecm_plan.data())
                                    / static_cast<double>(ecm_plan.size())),
                     stage_name);

            auto ecm_res = ecm_parallel(n, stage.curves, stage.B1, stage.B2, status_cb);
            if (ecm_res.has_value()) {
                T f = *ecm_res;
                auto sub1 = factor_parallel_generic(f, status_cb, profiler, sieve_method, progress_cb);
                auto sub2 = factor_parallel_generic(::bignum::ops::div_q(n, f),
                                                    status_cb,
                                                    profiler,
                                                    sieve_method,
                                                    progress_cb);
                factors.insert(factors.end(), sub1.begin(), sub1.end());
                factors.insert(factors.end(), sub2.begin(), sub2.end());
                return factors;
            }
        }
    }

    // 7. Решето по явному выбору пользователя — без ограничения по разрядности
    //    (выше kSieveHardMaxBits метод неприменим, это проверяется в try_sieve).
    if (sieve_method == kSieveForced && !::u128::Globals::LoadStop()) {
        if (try_sieve(true))
            return factors;
        if (status_cb)
            status_cb("SIQS не нашёл делитель: выбран принудительный SIQS, "
                      "другие решёта не запускаются.");
        if (profiler) profiler->finish();
        factors.push_back(n);
        return factors;
    }

    // 8. Решето в авторежиме — до kSieveMaxBits: выше время растёт слишком быстро.
    const bool sieve_allowed = (sieve_method == kSieveAuto && n_bits <= kSieveMaxBits);
    if (sieve_allowed && !::u128::Globals::LoadStop()) {
        if (try_sieve(false))
            return factors;
    }

    if (profiler) profiler->finish();
    factors.push_back(n);
    return factors;
}

// --- Генерация полупростого числа p×q (semiprime, как в криптографии) ---

namespace semiprime_detail {

// Простые числа < 5000 для дешёвого пробного деления кандидатов
inline const std::vector<uint32_t> &small_primes()
{
    static const std::vector<uint32_t> primes = [] {
        std::vector<bool> sieve(5000, true);
        sieve[0] = sieve[1] = false;
        for (uint32_t i = 2; i * i < 5000; ++i) {
            if (sieve[i]) {
                for (uint32_t j = i * i; j < 5000; j += i)
                    sieve[j] = false;
            }
        }
        std::vector<uint32_t> res;
        for (uint32_t i = 2; i < 5000; ++i) {
            if (sieve[i])
                res.push_back(i);
        }
        return res;
    }();
    return primes;
}

} // namespace semiprime_detail

// Генерирует полупростое число p×q ровно bitWidth бит: оба сомножителя простые,
// p длиной ровно bitWidth/2 бит, q — ровно bitWidth - bitWidth/2 бит (при нечётной
// bitWidth на один бит длиннее p), p != q.
// При срабатывании остановки (Globals::LoadStop) возвращает T(0).
// progress-колбэк получает число проверенных кандидатов (вызывается нечасто).
template <typename T>
T generate_semiprime(unsigned bitWidth,
                     T *p_out = nullptr,
                     T *q_out = nullptr,
                     const std::function<void(uint64_t)> &progress = nullptr,
                     uint64_t seed = 0)
{
    // Длины сомножителей. Для чётной bitWidth оба равны bitWidth/2. Для нечётной q
    // берётся на 1 бит длиннее p: при равных длинах произведение имеет максимум
    // bitWidth-1 бит, поэтому условие (p*q).bit_width() == bitWidth никогда не
    // выполнялось бы и поиск зацикливался (зависание генератора).
    const unsigned p_bits = bitWidth / 2;
    const unsigned q_bits = bitWidth - p_bits;

    // Защита от заведомо невыполнимых запросов: крошечная битность (в т.ч. нулевая,
    // из-за чего bits-1 переполнялся бы) либо битность, которая вместе с произведением
    // p*q не помещается в тип T (переполнение делало бы условие bit_width невыполнимым),
    // привели бы к вечному поиску — сразу возвращаем 0.
    if (bitWidth < 4 || p_bits == 0 || bitWidth > sizeof(T) * 8)
        return T(0);

    static thread_local std::mt19937_64 rng(std::random_device{}());
    if (seed != 0) // Детерминированный режим (тесты/стенд): воспроизводимое число
        rng.seed(seed);
    std::uniform_int_distribution<uint64_t> dist;

    // Случайное нечётное число с установленным старшим битом (ровно `bits` бит)
    auto random_odd = [&](unsigned bits) {
        T val(0);
        for (unsigned made = 0; made < bits; made += 64) {
            unsigned take = std::min(64u, bits - made);
            uint64_t chunk = dist(rng);
            if (take < 64)
                chunk &= (1ULL << take) - 1;

            val = (val << take) | T(chunk);
        }
        return val | (T(1) << (bits - 1)) | T(1);
    };

    // Дешёвый предварительный отсев: пробное деление по простым < 5000
    auto cheap_composite_check = [](const T &cand) {
        for (uint32_t sp : semiprime_detail::small_primes()) {
            if (cand <= T(sp))
                return false; // малое число — дальше решает Миллер—Рабин

            if (::bignum::ops::rem_u32(cand, sp) == 0)
                return true; // делится — составное
        }
        return false;
    };

    // Поиск простого длиной ровно `bits` бит
    auto find_prime = [&](T &out, unsigned bits) {
        for (;;) {
            if (::u128::Globals::LoadStop())
                return false;

            out = random_odd(bits);
            if (cheap_composite_check(out) || !is_prime_generic(out, 3))
                continue;

            return true;
        }
    };

    uint64_t candidates = 0;
    // Если у слишком маленькой p нет ни одного q с произведением нужной длины,
    // перебор завис бы на этой p навсегда — поэтому после такого числа неудачных
    // попыток p перерисовывается (гарантия выхода и отсутствие «зависших» хвостов).
    constexpr uint64_t kMaxTrialsPerP = 256;
    for (;;) {
        T p, q;
        if (!find_prime(p, p_bits))
            return T(0);

        // Ищем второе простое: q != p и p*q имеет ровно bitWidth бит
        bool found = false;
        for (uint64_t trial = 0; trial < kMaxTrialsPerP; ++trial) {
            if (::u128::Globals::LoadStop())
                return T(0);

            T cand;
            if (!find_prime(cand, q_bits))
                return T(0);

            ++candidates;
            if (progress && (candidates % 512) == 0)
                progress(candidates);

            if (cand == p)
                continue;

            if ((p * cand).bit_width() != bitWidth)
                continue; // p слишком мала для этой длины q — берём новую пару

            q = cand;
            found = true;
            break;
        }
        if (!found)
            continue;

        // Финальная проверка с высокой надёжностью (k=25, ложноположительное ~1e-15)
        if (!is_prime_generic(p, 25) || !is_prime_generic(q, 25))
            continue; // редчайшее ложное срабатывание фильтра — новая попытка

        if (progress)
            progress(candidates);

        if (p_out)
            *p_out = p;
        if (q_out)
            *q_out = q;

        return p * q;
    }
}

} // namespace bignum::generic_ecm
