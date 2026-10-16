// Частота структурной патологии (E_max>=2 с недостижимым ядром) по семействам
// матриц: считаем попытки до успеха и разбивку отказов решателя на Â.
#include "gnfs.hpp"
#include <iostream>
#include <random>
#include <algorithm>

using namespace bignum::gnfs::detail;

static Gf2Matrix make_planted_matrix(size_t M, size_t N, size_t nnz_row,
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
    for (size_t r = M - planted; r < M; ++r) {
        std::vector<uint8_t> acc(N, 0);
        size_t take = 2 + rng() % 5;
        for (size_t k = 0; k < take; ++k) {
            size_t src = rng() % r;
            for (uint32_t c : rows[src])
                acc[c] ^= 1;
        }
        for (uint32_t c = 0; c < N; ++c)
            if (acc[c])
                rows[r].push_back(c);
        if (rows[r].empty())
            rows[r].push_back(static_cast<uint32_t>(rng() % N));
    }
    for (size_t r = 0; r < M; ++r) {
        B.row_ptr[r] = static_cast<uint32_t>(B.row_cols.size());
        B.row_cols.insert(B.row_cols.end(), rows[r].begin(), rows[r].end());
    }
    B.row_ptr[M] = static_cast<uint32_t>(B.row_cols.size());
    B.build_csc();
    return B;
}

int main(int argc, char **argv)
{
    const size_t runs = argc > 1 ? std::stoul(argv[1]) : 20;
    const int max_attempts = argc > 2 ? std::stoi(argv[2]) : 300;
    const size_t cases[][4] = {
        // M, N, nnz_row, planted
        {120, 100, 6, 8},
        {120, 100, 12, 8},
        {120, 100, 40, 8},
        {120, 120, 40, 8},
        {400, 320, 40, 20},
        {400, 320, 8, 20},
    };
    for (const auto &cs : cases) {
        const size_t M = cs[0], N = cs[1], nnz_row = cs[2], planted = cs[3];
        std::mt19937_64 rng(0xBEEF + M * 131 + N * 17 + nnz_row);
        int pathological = 0;
        for (size_t run = 0; run < runs; ++run) {
            Gf2Matrix B = make_planted_matrix(M, N, nnz_row, planted, rng);
            int found = static_cast<int>(gf2_find_kernel_vectors(
                B, static_cast<size_t>(max_attempts), []() { return false; },
                [](const std::vector<uint64_t> &x) { return !x.empty(); }));
            if (!found)
                ++pathological;
        }
        std::cout << "M=" << M << " N=" << N << " nnz=" << nnz_row
                  << " planted=" << planted << ": патологичных " << pathological
                  << "/" << runs << " (лимит " << max_attempts << " попыток)\n";
    }
    return 0;
}
