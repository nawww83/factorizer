#pragma once

// ---------------------------------------------------------------------------
// siqs/relation.hpp — соотношения и частичные суммы по простым.
//
// FactorExp — пара (индекс простого, показатель); RawRelation — то, что
// выдаёт сито (большие простые, X, знак и разложение по базе);
// FullRelation — строка матрицы после склейки 1LP/2LP.
// ---------------------------------------------------------------------------

#include "common.hpp"

namespace bignum::siqs_xor {

// --- Соотношения -----------------------------------------------------------

struct FactorExp
{
    u32 idx;
    u32 cnt;
};

struct RawRelation
{
    // Большие простые, которых нет в базе (по возрастанию). Полное
    // соотношение — n_lp = 0; одно большое простое — n_lp = 1 и lp[0];
    // два — n_lp = 2 и lp[0], lp[1]. Третья ячейка оставлена под 3LP.
    std::array<u64, 3> lp{1, 1, 1};
    uint8_t n_lp = 0;
    Big X_val;
    int sign;
    std::vector<FactorExp> exps;
};

struct FullRelation
{
    Big X_val;
    int sign;
    std::vector<FactorExp> exps;
};

inline void merge_exps(std::vector<FactorExp> &dst, const std::vector<FactorExp> &src)
{
    for (const auto &fe : src)
    {
        bool found = false;
        for (auto &d : dst)
        {
            if (d.idx == fe.idx)
            {
                d.cnt += fe.cnt;
                found = true;
                break;
            }
        }
        if (!found)
            dst.push_back(fe);
    }
}

// То же для одного множителя: в горячем цикле разложения вызов через
// initializer_list создавал вектор (аллокация) на каждый найденный простой.
inline void merge_exp(std::vector<FactorExp> &dst, u32 idx, u32 cnt)
{
    for (auto &d : dst)
    {
        if (d.idx == idx)
        {
            d.cnt += cnt;
            return;
        }
    }
    dst.push_back({idx, cnt});
}

} // namespace bignum::siqs_xor
