// Тесты разреженного решателя СЛАУ GF(2) (метод Видеманна) из gnfs.hpp:
//   1) юнит-тест gf2_find_kernel_vectors на синтетических разреженных матрицах
//      с посаженными зависимостями (сверка B^T x == 0, x != 0);
//   2) сквозной прогон number_field_sieve на 128-битном полупростом числе.
//
// Сборка: g++ -std=c++20 -O2 -pthread -o gf2_la_test tests/gf2_la_test.cpp
// Запуск: ./gf2_la_test [только_юнит|только_e2e]

#include "../gnfs.hpp"

#include <iostream>
#include <random>

using bignum::gnfs::detail::Gf2Matrix;
using bignum::gnfs::detail::gf2_find_kernel_vectors;
using bignum::gnfs::detail::gf2_words;

namespace {

// Эталон: число нулевых строк после приведения Гауссом (левое ядро B^T).
size_t ref_left_nullity(const Gf2Matrix &B)
{
    const size_t M = B.rows_n;
    const size_t N = B.cols_n;
    const size_t W = (N + 63) / 64;
    std::vector<std::vector<uint64_t>> m(M, std::vector<uint64_t>(W, 0));
    for (size_t r = 0; r < M; ++r)
        for (uint32_t k = B.row_ptr[r]; k < B.row_ptr[r + 1]; ++k) {
            const uint32_t c = B.row_cols[k];
            m[r][c / 64] |= 1ULL << (c % 64);
        }
    size_t row = 0;
    for (size_t col = 0; col < N && row < M; ++col) {
        const uint64_t mask = 1ULL << (col % 64);
        size_t piv = row;
        while (piv < M && !(m[piv][col / 64] & mask))
            ++piv;
        if (piv == M)
            continue;
        std::swap(m[row], m[piv]);
        for (size_t r = 0; r < M; ++r)
            if (r != row && (m[r][col / 64] & mask))
                for (size_t w = col / 64; w < W; ++w)
                    m[r][w] ^= m[row][w];
        ++row;
    }
    size_t nullity = 0;
    for (size_t r = 0; r < M; ++r) {
        bool zero = true;
        for (size_t w = 0; w < W && zero; ++w)
            if (m[r][w])
                zero = false;
        nullity += zero ? 1 : 0;
    }
    return nullity;
}

// Проверка: B^T x == 0 и x != 0
bool check_dependency(const Gf2Matrix &B, const std::vector<uint64_t> &x)
{
    const size_t mw = gf2_words(B.rows_n);
    if (x.size() != mw)
        return false;
    bool nonzero = false;
    for (size_t i = 0; i < mw; ++i)
        nonzero = nonzero || (x[i] != 0);
    if (!nonzero)
        return false;
    std::vector<uint64_t> u(gf2_words(B.cols_n), 0);
    bignum::gnfs::detail::gf2_mul_Bt(B, x.data(), u.data());
    for (uint64_t w : u)
        if (w)
            return false;
    return true;
}

// Синтетическая матрица: M0 случайных разреженных строк + planted строки,
// являющиеся XOR'ом случайных подмножеств (гарантированные зависимости).
Gf2Matrix make_planted_matrix(size_t M, size_t N, size_t nnz_row,
                              size_t planted, std::mt19937_64 &rng)
{
    Gf2Matrix B;
    B.rows_n = M;
    B.cols_n = N;
    B.row_ptr.resize(M + 1);
    std::vector<std::vector<uint32_t>> rows(M);
    std::uniform_int_distribution<uint32_t> col_d(0, static_cast<uint32_t>(N - 1));
    for (size_t r = 0; r < M - planted; ++r) {
        std::vector<uint32_t> s;
        while (s.size() < nnz_row) {
            uint32_t c = col_d(rng);
            if (std::find(s.begin(), s.end(), c) == s.end())
                s.push_back(c);
        }
        std::sort(s.begin(), s.end());
        rows[r] = std::move(s);
    }
    // planted строки = XOR случайных подмножеств уже построенных строк
    for (size_t r = M - planted; r < M; ++r) {
        std::vector<uint8_t> acc(N, 0);
        size_t take = 2 + rng() % 5;
        for (size_t k = 0; k < take; ++k) {
            size_t src = rng() % r; // любая строка выше (в т.ч. planted)
            for (uint32_t c : rows[src])
                acc[c] ^= 1;
        }
        for (uint32_t c = 0; c < N; ++c)
            if (acc[c])
                rows[r].push_back(c);
        if (rows[r].empty())
            rows[r].push_back(static_cast<uint32_t>(rng() % N)); // не пустая строка
    }
    for (size_t r = 0; r < M; ++r)
        B.row_ptr[r] = static_cast<uint32_t>(B.row_cols.size()),
        B.row_cols.insert(B.row_cols.end(), rows[r].begin(), rows[r].end());
    B.row_ptr[M] = static_cast<uint32_t>(B.row_cols.size());
    B.build_csc();
    return B;
}

int unit_test()
{
    std::mt19937_64 rng(12345);
    const struct
    {
        size_t M, N, nnz, planted, runs;
    } cases[] = {
        {120, 100, 6, 8, 30},
        {300, 240, 8, 12, 20},
        {700, 500, 10, 20, 8},
    };
    int failures = 0;

    // Проверяем два слоя: (1) быстрый итерационный Видеманд — его доля успехов
    // на «типовых» матрицах должна быть высокой (регистресс-контроль);
    // (2) производственный контракт «Видеманд + точный fallback» — 100%.
    auto try_one = [](const Gf2Matrix &B, bool with_fallback) -> bool {
        std::vector<std::vector<uint64_t>> got;
        auto on_cand = [&](const std::vector<uint64_t> &x) {
            if (!check_dependency(B, x)) {
                std::cerr << "  [X] неверная зависимость B^T x = 0!\n";
                return true; // стоп — фиксируем провал
            }
            got.push_back(x);
            return got.size() >= 2; // пары независимых векторов достаточно
        };
        size_t found = gf2_find_kernel_vectors(B, 64, []() { return false; }, on_cand);
        if (!found && with_fallback)
            found = gf2_dense_kernel_combinations(B, 48, []() { return false; }, on_cand);
        return found >= 1 && got.size() >= 1 &&
               std::all_of(got.begin(), got.end(),
                           [&](const auto &x) { return check_dependency(B, x); });
    };

    for (const auto &cs : cases) {
        size_t ok_w = 0, ok_full = 0, tested = 0;
        size_t expect = 0;
        for (size_t run = 0; run < cs.runs; ++run) {
            Gf2Matrix B = make_planted_matrix(cs.M, cs.N, cs.nnz, cs.planted, rng);
            expect = ref_left_nullity(B); // должен быть >= planted
            if (expect == 0) {
                std::cerr << "  [!] матрица " << cs.M << "x" << cs.N
                          << ": эталонное ядро пустое — пропуск\n";
                continue;
            }
            ++tested;
            const bool w = try_one(B, false);
            ok_w += w ? 1 : 0;
            // полный контракт: Видеманд, при провале — точный плотный Гаусс
            ok_full += (w || try_one(B, true)) ? 1 : 0;
        }
        const bool w_pass = tested > 0 && ok_w * 5 >= tested * 3; // >=60% Видемандом
        const bool full_pass = (ok_full == tested);
        failures += (w_pass && full_pass) ? 0 : 1;
        std::cout << (w_pass && full_pass ? "  OK   " : "  FAIL ")
                  << "M=" << cs.M << " N=" << cs.N << " nnz=" << cs.nnz
                  << " planted=" << cs.planted << " (|left ker|~" << expect << "): "
                  << "Видеманд " << ok_w << "/" << tested
                  << ", с fallback " << ok_full << "/" << tested << "\n";
    }
    std::cout << (failures ? "ЮНИТ-ТЕСТ: ПРОВАЛ (" + std::to_string(failures) + ")\n"
                           : "ЮНИТ-ТЕСТ: OK\n");
    return failures;
}

int e2e_test()
{
    using T = ::bignum::U128;
    // 128-битное полупростое: (2^64-59) * (2^64-83)
    const T p("18446744073709551557");
    const T q("18446744073709551533");
    const T n = p * q;
    std::cout << "E2E: GNFS, n = " << n.toString() << " (" << n.bit_width() << " бит)\n";

    auto t0 = std::chrono::steady_clock::now();
    auto res = ::bignum::gnfs::number_field_sieve(n, [](const std::string &s) {
        std::cout << "  " << s << std::endl;
    });
    auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (!res.has_value()) {
        std::cout << "E2E: делитель не найден за " << sec << " c\n";
        return 1;
    }
    const T &g = *res;
    const bool good = (g > T(1)) && (g < n) && (n % g == T(0));
    std::cout << "E2E: делитель " << g.toString() << (good ? " — ВЕРНЫЙ" : " — ОШИБКА")
              << " (" << sec << " c)\n";
    return good ? 0 : 1;
}

} // namespace

int main(int argc, char **argv)
{
    const std::string mode = (argc > 1) ? argv[1] : "";
    int rc = 0;
    if (mode != "e2e")
        rc |= unit_test();
    if (mode != "unit")
        rc |= e2e_test();
    return rc;
}
