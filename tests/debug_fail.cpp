// Причины отказов попыток решателя на Â-подходе (тот же генератор матриц,
// что в юнит-тесте): считаем попытку до первого успеха на каждой из 30 матриц.
#include "gnfs.hpp"
#include <iostream>
#include <random>

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

int main()
{
    std::mt19937_64 rng(12345);
    const size_t M = 120, N = 100, nnz_row = 6, planted = 8;
    const size_t mw = gf2_words(M), nw = gf2_words(N);
    const size_t dw = mw + nw;
    const size_t padded = dw * 64;
    const size_t seq_len = 2 * padded + 2;

    int worst = 0;
    for (int run = 0; run < 30; ++run) {
        Gf2Matrix B = make_planted_matrix(M, N, nnz_row, planted, rng);

        if (run == 14 || run == 18 || run == 23) {
            std::cerr << "MATRIX " << run << " " << M << " " << N << "\n";
            for (size_t r = 0; r < M; ++r) {
                for (uint32_t k = B.row_ptr[r]; k < B.row_ptr[r + 1]; ++k)
                    std::cerr << B.row_cols[k] << " ";
                std::cerr << ";";
            }
            std::cerr << "\n";
        }
        Gf2Rng grng(0x1234 + run);
        int c_e0 = 0, c_x0 = 0, c_inv = 0, attempt_done = 0;
        std::vector<uint64_t> s(dw), t(dw), cur(dw), nxt(dw);
        std::vector<uint64_t> x(dw), acc(dw), scratch(nw), a(mw);
        std::vector<uint8_t> seq;
        seq.reserve(seq_len);

        for (int attempt = 1; attempt <= 500; ++attempt) {
            auto rand_vec = [&](std::vector<uint64_t> &v) {
                for (size_t i = 0; i < dw; ++i)
                    v[i] = grng.next();
                if (!gf2_any(v.data(), dw))
                    v[0] = 1;
            };
            rand_vec(s);
            rand_vec(t);
            seq.clear();
            std::copy(t.begin(), t.end(), cur.begin());
            for (size_t k = 0; k < seq_len; ++k) {
                seq.push_back(gf2_dot(s.data(), cur.data(), dw) ? 1 : 0);
                gf2_mul_A(B, cur.data(), nxt.data());
                cur.swap(nxt);
            }
            auto mu = berlekamp_massey(seq);
            size_t e = 0;
            while (e < mu.size() && mu[e] == 0)
                ++e;
            if (e == 0) {
                ++c_e0;
                continue;
            }
            std::fill(acc.begin(), acc.end(), 0);
            const size_t deg_h = mu.size() - 1 - e;
            for (size_t i = deg_h + 1; i-- > 0;) {
                if (gf2_any(acc.data(), dw)) {
                    gf2_mul_A(B, acc.data(), nxt.data());
                    acc.swap(nxt);
                }
                if (mu[e + i])
                    gf2_xor(acc.data(), t.data(), dw);
            }
            x = acc;
            for (size_t i = 1; i < e; ++i) {
                gf2_mul_A(B, x.data(), nxt.data());
                x.swap(nxt);
            }
            std::copy(x.begin(), x.begin() + mw, a.begin());
            const size_t tail = M & 63;
            if (tail)
                a[mw - 1] &= (1ULL << tail) - 1;
            if (!gf2_any(a.data(), mw)) {
                ++c_x0;
                continue;
            }
            gf2_mul_Bt(B, a.data(), scratch.data());
            if (gf2_any(scratch.data(), nw)) {
                ++c_inv; // для Â «невозможный» случай — контроль теории
                continue;
            }
            attempt_done = attempt;
            break;
        }
        worst = std::max(worst, attempt_done);
        std::cout << "run " << run << ": попыток=" << attempt_done
                  << " (e0=" << c_e0 << " x0=" << c_x0 << " inv=" << c_inv
                  << ")" << (attempt_done ? "" : "   <-- ПРОВАЛ") << "\n";
    }
    std::cout << "худшее число попыток: " << worst << "\n";
    return 0;
}
