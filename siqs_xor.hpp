#pragma once

// ---------------------------------------------------------------------------
// siqs_xor.hpp — самодостаточное квадратичное решето (SIQS) с 2LP-циклами.
//
// Портировано из tests/siqs_xor.cpp: одиночный проход по A-полиномам,
// L1/L2-сито по байтовым логарифмам, hit-map для больших простых, лес
// циклов 2LP и разреженный гаусс над GF(2) с XOR-историей строк.
//
// Отличия от эталонного файла:
//   * нет OpenMP — батчи A-полиномов считаются std::thread-воркерами
//     (используется cpu_topology, как в ECM);
//   * весь вывод идёт через status_cb (std::cout в приложении не нужен);
//   * участвует глобальный флаг остановки u128::Globals::LoadStop();
//   * вместо собственного U512 используется фиксированный bignum::UInt<512>.
// ---------------------------------------------------------------------------

#include "bignum.hpp"
#include "cpu_topology.hpp"
#include "u128_utils.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(__SSE2__)
#  include <emmintrin.h>
#endif

namespace bignum::siqs_xor {

using u32 = uint32_t;
using u64 = uint64_t;
using i64 = int64_t;
using Big = UInt<512>;

inline bool stop_requested()
{
    return ::u128::Globals::LoadStop();
}

// --- 64-битная модулярная арифметика ---------------------------------------

inline u64 mul_mod64(u64 a, u64 b, u64 m)
{
#if defined(__SIZEOF_INT128__)
    return static_cast<u64>((static_cast<unsigned __int128>(a) * b) % m);
#else
    u64 hi, lo;
    intrinsics::mul64_to_128(a, b, hi, lo);
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
    u64 rem = 0;
    _udiv128(hi, lo, m, &rem);
    return rem;
#else
    u64 rem = 0;
    for (int bit = 127; bit >= 0; --bit)
    {
        const u64 bit_val = (bit >= 64) ? ((hi >> (bit - 64)) & 1ULL) : ((lo >> bit) & 1ULL);
        const u64 two = (rem << 1) | bit_val;
        if ((rem >> 63) || two >= m)
            two -= m;
        rem = two;
    }
    return rem;
#endif
#endif
}

inline u64 pow_mod64(u64 base, u64 exp, u64 m)
{
    if (m == 1)
        return 0;
    u64 res = 1;
    base %= m;
    while (exp > 0)
    {
        if (exp & 1)
            res = mul_mod64(res, base, m);
        base = mul_mod64(base, base, m);
        exp >>= 1;
    }
    return res;
}

inline u64 gcd64(u64 a, u64 b)
{
    while (b)
    {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

inline u64 inv_mod64(u64 a, u64 m)
{
    i64 t = 0, newt = 1;
    i64 r = static_cast<i64>(m), newr = static_cast<i64>(a % m);
    while (newr != 0)
    {
        i64 q = r / newr;
        i64 tmp_t = t - q * newt;
        t = newt;
        newt = tmp_t;
        i64 tmp_r = r - q * newr;
        r = newr;
        newr = tmp_r;
    }
    if (t < 0)
        t += static_cast<i64>(m);
    return static_cast<u64>(t);
}

// Детерминированный тест Миллера — Рабина для u64.
inline bool is_prime64(u64 n)
{
    if (n < 2)
        return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL,
                  23ULL, 29ULL, 31ULL, 37ULL})
    {
        if (n == p)
            return true;
        if (n % p == 0)
            return false;
    }
    u64 d = n - 1;
    const int r = static_cast<int>(std::countr_zero(d));
    d >>= r;
    for (u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL,
                  23ULL, 29ULL, 31ULL, 37ULL})
    {
        if (n <= a)
            break;
        u64 x = pow_mod64(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool composite = true;
        for (int i = 0; i < r - 1; ++i)
        {
            x = mul_mod64(x, x, n);
            if (x == n - 1)
            {
                composite = false;
                break;
            }
        }
        if (composite)
            return false;
    }
    return true;
}

// Метод Брента — Полларда для 64-битного кофактора.
inline bool brent_factor_two_primes(u64 n, u64 &f1, u64 &f2)
{
    if ((n & 1) == 0)
    {
        f1 = 2;
        f2 = n >> 1;
        return true;
    }
    for (u64 c : {1ULL, 3ULL, 5ULL})
    {
        u64 y = 2, m = 32, g = 1, r = 1, q = 1, ys = 2, x = 2;
        while (g == 1)
        {
            x = y;
            for (u64 i = 0; i < r; ++i)
                y = (mul_mod64(y, y, n) + c) % n;
            u64 k = 0;
            while (k < r && g == 1)
            {
                ys = y;
                const u64 lim = std::min(m, r - k);
                for (u64 i = 0; i < lim; ++i)
                {
                    y = (mul_mod64(y, y, n) + c) % n;
                    const u64 diff = (x > y) ? (x - y) : (y - x);
                    q = mul_mod64(q, diff, n);
                }
                g = gcd64(q, n);
                k += m;
            }
            r <<= 1;
            if (r > 1024)
                break;
        }
        if (g == n)
        {
            while (true)
            {
                ys = (mul_mod64(ys, ys, n) + c) % n;
                const u64 diff = (x > ys) ? (x - ys) : (ys - x);
                g = gcd64(diff, n);
                if (g > 1)
                    break;
            }
        }
        if (g > 1 && g < n)
        {
            f1 = g;
            f2 = n / g;
            return true;
        }
    }
    return false;
}

// Алгоритм Тонелли — Шэнкса: r^2 == n (mod p).
inline u32 tonelli_shanks(u64 n, u32 p)
{
    n %= p;
    if (n == 0)
        return 0;
    if (p == 2)
        return static_cast<u32>(n);
    if ((p & 3) == 3)
        return static_cast<u32>(pow_mod64(n, (p + 1) / 4, p));
    u32 q = p - 1, s = 0;
    while ((q & 1) == 0)
    {
        q >>= 1;
        ++s;
    }
    u32 z = 2;
    while (pow_mod64(z, (p - 1) / 2, p) != p - 1)
        ++z;
    u32 m = s;
    u64 c = pow_mod64(z, q, p);
    u64 t = pow_mod64(n, q, p);
    u64 r = pow_mod64(n, (q + 1) / 2, p);
    while (t != 1)
    {
        u32 i = 1;
        u64 temp = (t * t) % p;
        while (temp != 1)
        {
            temp = (temp * temp) % p;
            ++i;
        }
        const u64 b = pow_mod64(c, 1ULL << (m - i - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }
    return static_cast<u32>(r);
}

// --- Модулярные операции над 512-битным типом ------------------------------

inline Big pow_mod_big(Big base, u64 exp, const Big &m)
{
    BarrettReducer<512> br(m);
    Big res(1);
    if (base >= m)
        base = base % m;
    while (exp > 0)
    {
        if (exp & 1)
            res = br.mul(res, base);
        exp >>= 1;
        if (exp != 0)
            base = br.sqr(base);
    }
    return res;
}

inline Big inv_mod_big(Big a, const Big &m)
{
    a = a % m;
    Big u = a, v = m;
    Big x1(1), x2(0);
    while (!u.is_zero())
    {
        Big q, r;
        q = Big::divide(v, u, &r);
        v = u;
        u = r;
        Big qx1 = Big::mul_mod(q, x1, m);
        Big nx = (x2 >= qx1) ? (x2 - qx1) : (m - (qx1 - x2));
        x2 = x1;
        x1 = nx;
    }
    return x2;
}

// --- Соотношения -----------------------------------------------------------

struct FactorExp
{
    u32 idx;
    u32 cnt;
};

struct RawRelation
{
    u64 u_lp;
    u64 v_lp;
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

// Граф циклов 2LP на ориентированном остовном лесе.
class LargePrimeGraph
{
    struct TreeEdge
    {
        u64 parent_node;
        Big X_val;
        int sign;
        std::vector<FactorExp> exps;
        u64 p1, p2;
    };

    Big M;
    std::unordered_map<u64, u64> uf_parent;
    std::unordered_map<u64, TreeEdge> tree_parent;

    u64 find_root(u64 v)
    {
        u64 curr = v;
        while (true)
        {
            auto it = uf_parent.find(curr);
            if (it == uf_parent.end() || it->second == curr)
                break;
            curr = it->second;
        }
        const u64 root = curr;
        curr = v;
        while (true)
        {
            auto it = uf_parent.find(curr);
            if (it == uf_parent.end() || it->second == root)
                break;
            const u64 nxt = it->second;
            it->second = root;
            curr = nxt;
        }
        if (uf_parent.find(root) == uf_parent.end())
            uf_parent[root] = root;
        return root;
    }

    void reroot_tree(u64 node)
    {
        std::vector<std::pair<u64, TreeEdge>> path;
        u64 curr = node;
        while (true)
        {
            auto it = tree_parent.find(curr);
            if (it == tree_parent.end())
                break;
            path.push_back({curr, it->second});
            curr = it->second.parent_node;
        }
        for (int i = static_cast<int>(path.size()) - 1; i >= 0; --i)
        {
            const u64 child = path[i].first;
            TreeEdge te = path[i].second;
            const u64 par = te.parent_node;
            te.parent_node = child;
            tree_parent[par] = std::move(te);
        }
        tree_parent.erase(node);
    }

public:
    u64 num_edges = 0;
    u64 cycles_found = 0;

    explicit LargePrimeGraph(const Big &m_val) : M(m_val) {}

    bool add_edge_and_extract_cycle(const RawRelation &rel, FullRelation &out_rel)
    {
        const u64 u = rel.u_lp, v = rel.v_lp;
        if (u == 1 && v == 1)
        {
            out_rel = {rel.X_val, rel.sign, rel.exps};
            return true;
        }
        if (u == v)
        {
            const Big inv_u = inv_mod_big(Big(u), M);
            out_rel = {Big::mul_mod(rel.X_val, inv_u, M), rel.sign, rel.exps};
            return true;
        }

        const u64 ru = find_root(u);
        const u64 rv = find_root(v);
        const u64 ep1 = (u > 1) ? u : 0;
        const u64 ep2 = (v > 1) ? v : 0;

        if (ru != rv)
        {
            reroot_tree(u);
            tree_parent[u] = {v, rel.X_val, rel.sign, rel.exps, ep1, ep2};
            uf_parent[ru] = rv;
            ++num_edges;
            return false;
        }

        std::unordered_set<u64> anc_u;
        u64 curr = u;
        while (true)
        {
            anc_u.insert(curr);
            auto it = tree_parent.find(curr);
            if (it == tree_parent.end())
                break;
            curr = it->second.parent_node;
        }
        curr = v;
        while (anc_u.find(curr) == anc_u.end())
            curr = tree_parent[curr].parent_node;
        const u64 lca = curr;

        Big comb_X = rel.X_val;
        int comb_sign = rel.sign;
        std::vector<FactorExp> comb_exps = rel.exps;
        std::unordered_map<u64, u32> lp_counts;
        if (ep1)
            ++lp_counts[ep1];
        if (ep2)
            ++lp_counts[ep2];

        for (u64 start_node : {u, v})
        {
            curr = start_node;
            while (curr != lca)
            {
                const auto &te = tree_parent[curr];
                comb_X = Big::mul_mod(comb_X, te.X_val, M);
                comb_sign *= te.sign;
                merge_exps(comb_exps, te.exps);
                if (te.p1)
                    ++lp_counts[te.p1];
                if (te.p2)
                    ++lp_counts[te.p2];
                curr = te.parent_node;
            }
        }

        Big denom(1);
        for (const auto &kv : lp_counts)
        {
            const Big lp_pow = pow_mod_big(Big(kv.first), kv.second / 2, M);
            denom = Big::mul_mod(denom, lp_pow, M);
        }
        const Big inv_denom = inv_mod_big(denom, M);
        comb_X = Big::mul_mod(comb_X, inv_denom, M);

        ++cycles_found;
        out_rel = {comb_X, comb_sign, std::move(comb_exps)};
        return true;
    }
};

// Разреженный гаусс над GF(2) с XOR-историей строк.
class MarkowitzGaussianGF2
{
public:
    std::vector<FullRelation> relations;
    std::unordered_map<u32, u32> col_freq;
    u32 multi_hit_primes = 0;
    u32 singleton_primes = 0;
    size_t last_check_count = 0;
    std::function<void(const std::string &)> status_cb;
    // Прогресс исключения (0..1 внутри фазы Гаусса), вызывается из основного потока.
    std::function<void(double, const std::string &)> progress_cb;

    void add_relation(FullRelation &&rel)
    {
        for (const auto &fe : rel.exps)
        {
            if (fe.cnt & 1)
            {
                const u32 prev = col_freq[fe.idx]++;
                if (prev == 0)
                    ++singleton_primes;
                else if (prev == 1)
                {
                    --singleton_primes;
                    ++multi_hit_primes;
                }
            }
        }
        relations.push_back(std::move(rel));
    }

    bool ready_to_check_rank() const
    {
        const size_t n_rels = relations.size();
        if (n_rels < last_check_count || n_rels - last_check_count < 500)
            return false;
        const size_t total_active = static_cast<size_t>(multi_hit_primes) + singleton_primes;
        return n_rels + 120 >= total_active;
    }

    std::vector<std::vector<u32>> solve_pruned_matrix()
    {
        last_check_count = relations.size();
        struct ActiveRow
        {
            std::unordered_set<u32> cols;
            int s_bit;
            std::vector<u32> history;
        };

        std::unordered_map<u32, ActiveRow> active_rows;
        std::unordered_map<u32, std::unordered_set<u32>> col_to_rows;

        for (u32 r_idx = 0; r_idx < static_cast<u32>(relations.size()); ++r_idx)
        {
            std::unordered_set<u32> odd_cols;
            for (const auto &fe : relations[r_idx].exps)
            {
                if (fe.cnt & 1)
                {
                    odd_cols.insert(fe.idx);
                    col_to_rows[fe.idx].insert(r_idx);
                }
            }
            const int s_bit = (relations[r_idx].sign == -1) ? 1 : 0;
            active_rows[r_idx] = {std::move(odd_cols), s_bit, {r_idx}};
        }

        auto remove_row = [&](u32 r_del)
        {
            auto it = active_rows.find(r_del);
            if (it == active_rows.end())
                return;
            for (u32 c : it->second.cols)
            {
                auto cit = col_to_rows.find(c);
                if (cit != col_to_rows.end())
                {
                    cit->second.erase(r_del);
                    if (cit->second.empty())
                        col_to_rows.erase(cit);
                }
            }
            active_rows.erase(it);
        };

        auto xor_history = [](std::vector<u32> &dst, const std::vector<u32> &src)
        {
            std::vector<u32> res;
            res.reserve(dst.size() + src.size());
            size_t i = 0, j = 0;
            while (i < dst.size() && j < src.size())
            {
                if (dst[i] < src[j])
                    res.push_back(dst[i++]);
                else if (src[j] < dst[i])
                    res.push_back(src[j++]);
                else
                {
                    ++i;
                    ++j;
                }
            }
            while (i < dst.size())
                res.push_back(dst[i++]);
            while (j < src.size())
                res.push_back(src[j++]);
            dst = std::move(res);
        };

        auto xor_row_into = [&](u32 r_src, u32 r_dst)
        {
            const auto &src = active_rows[r_src];
            auto &dst = active_rows[r_dst];
            for (u32 c : src.cols)
            {
                if (dst.cols.count(c))
                {
                    dst.cols.erase(c);
                    col_to_rows[c].erase(r_dst);
                    if (col_to_rows[c].empty())
                        col_to_rows.erase(c);
                }
                else
                {
                    dst.cols.insert(c);
                    col_to_rows[c].insert(r_dst);
                }
            }
            dst.s_bit ^= src.s_bit;
            xor_history(dst.history, src.history);
        };

        auto prune_singletons = [&]()
        {
            std::vector<u32> q;
            for (const auto &kv : col_to_rows)
                if (kv.second.size() == 1)
                    q.push_back(kv.first);
            while (!q.empty())
            {
                const u32 c = q.back();
                q.pop_back();
                auto cit = col_to_rows.find(c);
                if (cit != col_to_rows.end() && cit->second.size() == 1)
                {
                    const u32 r_del = *cit->second.begin();
                    std::vector<u32> aff(active_rows[r_del].cols.begin(),
                                         active_rows[r_del].cols.end());
                    remove_row(r_del);
                    for (u32 ca : aff)
                    {
                        auto ait = col_to_rows.find(ca);
                        if (ait != col_to_rows.end() && ait->second.size() == 1)
                            q.push_back(ca);
                    }
                }
            }
        };

        prune_singletons();
        if (active_rows.size() <= col_to_rows.size() + 5)
            return {};
        const size_t orig_core = active_rows.size();

        // Сжатие SGE для весов 2..5.
        for (size_t max_w : {2, 2, 3, 3, 4, 5})
        {
            std::vector<u32> cands;
            for (const auto &kv : col_to_rows)
                if (kv.second.size() >= 2 && kv.second.size() <= max_w)
                    cands.push_back(kv.first);
            for (u32 c : cands)
            {
                auto cit = col_to_rows.find(c);
                if (cit == col_to_rows.end() || cit->second.size() < 2
                    || cit->second.size() > max_w)
                    continue;
                if (active_rows.size() <= col_to_rows.size() + 24)
                    break;
                std::vector<u32> r_list(cit->second.begin(), cit->second.end());
                std::sort(r_list.begin(), r_list.end(), [&](u32 a, u32 b)
                          { return active_rows[a].cols.size() < active_rows[b].cols.size(); });
                const u32 pivot_r = r_list[0];
                if (active_rows[pivot_r].cols.size() > 40)
                    continue;
                for (size_t i = 1; i < r_list.size(); ++i)
                    xor_row_into(pivot_r, r_list[i]);
                remove_row(pivot_r);
            }
            prune_singletons();
        }

        if (status_cb)
            status_cb("SIQS: Гаусс: ядро после синглетонов " + std::to_string(orig_core)
                      + " -> сжато SGE до " + std::to_string(active_rows.size())
                      + " строк и " + std::to_string(col_to_rows.size()) + " столбцов");

        std::vector<u32> sorted_cols;
        sorted_cols.reserve(col_to_rows.size());
        for (const auto &kv : col_to_rows)
            sorted_cols.push_back(kv.first);
        std::sort(sorted_cols.begin(), sorted_cols.end(), [&](u32 a, u32 b)
                  { return col_to_rows[a].size() > col_to_rows[b].size(); });
        std::unordered_map<u32, u32> col_map;
        for (u32 i = 0; i < static_cast<u32>(sorted_cols.size()); ++i)
            col_map[sorted_cols[i]] = i;
        const u32 sign_bit = static_cast<u32>(sorted_cols.size());
        const u32 num_bits = sign_bit + 1;
        const size_t vec_words = (num_bits + 63) / 64;

        std::vector<u32> row_order;
        row_order.reserve(active_rows.size());
        for (const auto &kv : active_rows)
            row_order.push_back(kv.first);
        std::sort(row_order.begin(), row_order.end(), [&](u32 a, u32 b)
                  { return active_rows[a].cols.size() < active_rows[b].cols.size(); });

        const size_t hist_words = (row_order.size() + 63) / 64;
        struct PivotEntry
        {
            std::vector<u64> vec;
            std::vector<u64> hist;
        };
        std::unordered_map<int, PivotEntry> pivot_map;
        std::vector<std::vector<u32>> deps;

        for (size_t ord_i = 0; ord_i < row_order.size(); ++ord_i)
        {
            if ((ord_i & 0x3FF) == 0 && stop_requested())
                return deps;

            if (progress_cb && (ord_i & 0x3FF) == 0 && !row_order.empty())
                progress_cb(static_cast<double>(ord_i) / static_cast<double>(row_order.size()),
                            "исключение Гаусса");

            const u32 r_idx = row_order[ord_i];
            const auto &ar = active_rows[r_idx];
            std::vector<u64> vec(vec_words, 0);
            for (u32 c : ar.cols)
            {
                auto mit = col_map.find(c);
                if (mit != col_map.end())
                {
                    const u32 bit = mit->second;
                    vec[bit >> 6] |= (1ULL << (bit & 63));
                }
            }
            if (ar.s_bit)
                vec[sign_bit >> 6] |= (1ULL << (sign_bit & 63));

            std::vector<u64> hist(hist_words, 0);
            hist[ord_i >> 6] |= (1ULL << (ord_i & 63));
            int msb = -1;
            for (int w = static_cast<int>(vec_words) - 1; w >= 0; --w)
            {
                if (vec[w] != 0)
                {
                    msb = w * 64 + (63 - static_cast<int>(std::countl_zero(vec[w])));
                    break;
                }
            }
            while (msb >= 0)
            {
                auto pit = pivot_map.find(msb);
                if (pit != pivot_map.end())
                {
                    const auto &pv = pit->second.vec;
                    const auto &ph = pit->second.hist;
                    const int top_w = msb >> 6;
                    for (int w = 0; w <= top_w; ++w)
                        vec[w] ^= pv[w];
                    for (size_t w = 0; w < hist_words; ++w)
                        hist[w] ^= ph[w];
                    msb = -1;
                    for (int w = top_w; w >= 0; --w)
                    {
                        if (vec[w] != 0)
                        {
                            msb = w * 64 + (63 - static_cast<int>(std::countl_zero(vec[w])));
                            break;
                        }
                    }
                }
                else
                {
                    pivot_map[msb] = {std::move(vec), std::move(hist)};
                    break;
                }
            }
            if (msb < 0)
            {
                std::vector<u32> comb;
                for (size_t w = 0; w < hist_words; ++w)
                {
                    u64 word = hist[w];
                    while (word)
                    {
                        const int b = static_cast<int>(std::countr_zero(word));
                        const size_t idx_in_ord = w * 64 + b;
                        xor_history(comb, active_rows[row_order[idx_in_ord]].history);
                        word &= (word - 1);
                    }
                }
                if (!comb.empty())
                {
                    deps.push_back(std::move(comb));
                    if (deps.size() >= 16)
                        break;
                }
            }
        }
        return deps;
    }
};

// --- Параметры и конфигурация ----------------------------------------------

struct SIQSConfig
{
    Big M;
    Big kM;
    std::vector<u32> prime_base;
    std::unordered_map<u64, u32> prime_to_idx;
    std::vector<u32> root1;
    std::vector<u32> root2;
    std::vector<uint8_t> log_primes;
    u32 num_primes = 0;
    u64 last_prime = 0;
    u64 min_2lp_comp = 0;
    u64 lp_bound = 0;
    u64 lp2_bound = 0;
    u32 B_size = 0;
    u32 half_B = 0;
    int s = 0;
    std::vector<u32> pool_indices;
    u32 small_prime_skip = 0;
    u32 split_idx = 0;
};

// Сдвиг корней при переходе к следующему полиному: x += delta (add = true) или
// x -= delta (add = false) по модулю p. Условие `if (inv_A_mod_p[idx])` убрано:
// у выбранных в A простых delta = 0, поэтому сдвиг нулевой в обоих режимах.
// Скалярная ветка на каждом простом стоила дороже самой арифметики.
inline void update_roots(u32 num_primes, const u32 *__restrict prime,
                         const u32 *__restrict delta, u32 *__restrict sol1,
                         u32 *__restrict sol2, bool add)
{
    u32 idx = 1;
#if defined(__SSE2__)
    if (add)
    {
        // sol < p и delta < p, значит сумма < 2p < 2^31: знаковое сравнение
        // корректно, вычитание p делается по маске.
        for (; idx + 4 <= num_primes; idx += 4)
        {
            const __m128i d = _mm_loadu_si128(reinterpret_cast<const __m128i *>(delta + idx));
            const __m128i pp = _mm_loadu_si128(reinterpret_cast<const __m128i *>(prime + idx));
            __m128i a = _mm_add_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol1 + idx)), d);
            __m128i b = _mm_add_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol2 + idx)), d);
            a = _mm_sub_epi32(a, _mm_andnot_si128(_mm_srai_epi32(_mm_sub_epi32(a, pp), 31), pp));
            b = _mm_sub_epi32(b, _mm_andnot_si128(_mm_srai_epi32(_mm_sub_epi32(b, pp), 31), pp));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol1 + idx), a);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol2 + idx), b);
        }
    }
    else
    {
        for (; idx + 4 <= num_primes; idx += 4)
        {
            const __m128i d = _mm_loadu_si128(reinterpret_cast<const __m128i *>(delta + idx));
            const __m128i pp = _mm_loadu_si128(reinterpret_cast<const __m128i *>(prime + idx));
            __m128i a = _mm_sub_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol1 + idx)), d);
            __m128i b = _mm_sub_epi32(
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(sol2 + idx)), d);
            a = _mm_add_epi32(a, _mm_and_si128(pp, _mm_srai_epi32(a, 31)));
            b = _mm_add_epi32(b, _mm_and_si128(pp, _mm_srai_epi32(b, 31)));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol1 + idx), a);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(sol2 + idx), b);
        }
    }
#endif
    if (add)
    {
        for (; idx < num_primes; ++idx)
        {
            const u32 p = prime[idx];
            const u32 d = delta[idx];
            const u32 s1_n = sol1[idx] + d;
            sol1[idx] = (s1_n >= p) ? (s1_n - p) : s1_n;
            const u32 s2_n = sol2[idx] + d;
            sol2[idx] = (s2_n >= p) ? (s2_n - p) : s2_n;
        }
    }
    else
    {
        for (; idx < num_primes; ++idx)
        {
            const u32 p = prime[idx];
            const u32 d = delta[idx];
            sol1[idx] = (sol1[idx] >= d) ? (sol1[idx] - d) : (sol1[idx] + p - d);
            sol2[idx] = (sol2[idx] >= d) ? (sol2[idx] - d) : (sol2[idx] + p - d);
        }
    }
}

// Один батч A-полиномов (полностью однопоточный).
inline std::vector<RawRelation> process_one_A_batch(const SIQSConfig &cfg, u64 seed,
                                                    int poly_limit = 0)
{
    const u32 num_primes = cfg.num_primes;
    const u32 B_size = cfg.B_size;
    const u32 half_B = cfg.half_B;
    const int s = cfg.s;

    std::mt19937_64 rng(seed);
    std::vector<u32> chosen_indices;
    std::vector<u32> pool_copy = cfg.pool_indices;
    for (int i = 0; i < s; ++i)
    {
        std::uniform_int_distribution<size_t> dist(i, pool_copy.size() - 1);
        const size_t pick = dist(rng);
        std::swap(pool_copy[i], pool_copy[pick]);
        chosen_indices.push_back(pool_copy[i]);
    }
    std::sort(chosen_indices.begin(), chosen_indices.end());

    std::vector<uint8_t> is_chosen(num_primes, 0);
    Big A(1);
    for (u32 idx : chosen_indices)
    {
        is_chosen[idx] = 1;
        A = A * Big(cfg.prime_base[idx]);
    }

    std::vector<Big> B_terms(s);
    for (int j = 0; j < s; ++j)
    {
        const u32 idx = chosen_indices[j];
        const u64 q_j = cfg.prime_base[idx];
        const u64 r_j = cfg.root1[idx];
        Big A_div_q = A;
        A_div_q.divrem_u64(q_j);
        const u64 rem_a = A_div_q.mod_u64(q_j);
        const u64 inv = inv_mod64(rem_a, q_j);
        u64 gamma = mul_mod64(r_j, inv, q_j);
        if (gamma > q_j / 2)
            gamma = q_j - gamma;
        B_terms[j] = A_div_q * Big(gamma);
    }

    // B = sum(+-B_terms[j]) — держим модуль и знак отдельно.
    Big B_mag(0);
    bool B_neg = false;
    for (int j = 0; j < s; ++j)
        B_mag.add(B_terms[j]);
    auto b_signed_add = [&](const Big &term, bool term_neg)
    {
        if (B_neg == term_neg)
        {
            B_mag = B_mag + term;
        }
        else if (B_mag >= term)
        {
            B_mag = B_mag - term;
        }
        else
        {
            B_mag = term - B_mag;
            B_neg = term_neg;
        }
        if (B_mag.is_zero())
            B_neg = false;
    };

    std::vector<u32> inv_A_mod_p(num_primes, 0);
    std::vector<std::vector<u32>> B_terms_mod_p(s, std::vector<u32>(num_primes, 0));
    std::vector<u32> sol1(num_primes, 0);
    std::vector<u32> sol2(num_primes, 0);

    // floor(2^32 / p) для быстрого остатка в пробном делении: на 300 битах
    // s_idx % p считается для ~6 тыс. простых на каждого кандидата, а деление
    // на переменный делитель — это 20+ тактов. Умножение с восстановлением
    // остатка стоит ~5 тактов и даёт тот же результат.
    std::vector<u32> p_recip(num_primes, 0);

    for (u32 idx = 1; idx < num_primes; ++idx)
    {
        const u32 p = cfg.prime_base[idx];
        p_recip[idx] = static_cast<u32>((1ULL << 32) / p);
        if (is_chosen[idx])
            continue;
        const u64 a_mod = A.mod_u64(p);
        const u32 inv_A = static_cast<u32>(inv_mod64(a_mod, p));
        inv_A_mod_p[idx] = inv_A;
        for (int j = 0; j < s; ++j)
        {
            const u64 b_mod = B_terms[j].mod_u64(p);
            B_terms_mod_p[j][idx] = static_cast<u32>((2ULL * b_mod * inv_A) % p);
        }

        u64 B_mod = B_mag.mod_u64(p);
        if (B_neg && B_mod != 0)
            B_mod = static_cast<u64>(p) - B_mod;
        const u64 r1 = cfg.root1[idx];
        const u64 r2 = cfg.root2[idx];
        sol1[idx] = static_cast<u32>((((p + r1 - B_mod) % p) * static_cast<u64>(inv_A)
                                      + static_cast<u64>(half_B)) % p);
        sol2[idx] = static_cast<u32>((((p + r2 - B_mod) % p) * static_cast<u64>(inv_A)
                                      + static_cast<u64>(half_B)) % p);
    }

    std::vector<int> signs(s, 1);
    const int polys_per_A = 1 << (s - 1);
    // Ограничение числа полиномов на один A (poly_limit) позволяет набирать
    // соотношения малыми порциями: на небольших числах один полный батч даёт
    // десятки тысяч соотношений, и разреженный Гаусс потом считается впустую.
    const int polys_to_run = (poly_limit <= 0 || poly_limit > polys_per_A)
                                 ? polys_per_A
                                 : poly_limit;
    std::vector<RawRelation> found_relations;

    std::vector<uint8_t> sieve_logs(B_size);
    // Попадания крупных простых в блок: старшие 32 бита — позиция, младшие —
    // индекс простого. Раньше здесь лежали индексы, и второй проход (отбор
    // выживших) заново читал sol1/sol2 по случайным адресам: на 300 битах это
    // десятки тысяч промахов кэша на полином. Теперь проход идёт по плотному
    // массиву пар, а индексация нужна только для редких выживших.
    std::vector<u64> in_block_large1;
    std::vector<u64> in_block_large2;
    in_block_large1.reserve(4096);
    in_block_large2.reserve(4096);

    // Битовая карта выживших позиций и список попаданий крупных простых.
    // Хеш-таблица на полином здесь была лишней: выживших единицы, а проверок —
    // десятки тысяч, и каждая стоила похода в unordered_map. Карта на блок
    // в 128 КБ занимает 16 КБ и целиком лежит в L1.
    std::vector<u64> survivor_bits((B_size + 63) / 64);
    std::vector<u32> survivor_idx;
    // Отобранные пары (позиция, индекс простого) вместо двух параллельных
    // векторов: одна загрузка на проверку в цикле разложения.
    std::vector<u64> hit_pairs;
    std::vector<u32> large_divs;
    survivor_idx.reserve(1024);
    hit_pairs.reserve(64);
    large_divs.reserve(16);
    // Переиспользуемый буфер множителей: раньше вектор создавался на каждого
    // выжившего кандидата (сотни тысяч аллокаций за проход).
    std::vector<FactorExp> exps;
    exps.reserve(64);

    const int base_target_bits = static_cast<int>(cfg.kM.bit_width() / 2
                                                  + std::countr_zero(half_B));
    const int lp2_bits = static_cast<int>(std::bit_width(cfg.lp2_bound));
    const uint8_t target_log = static_cast<uint8_t>(
        std::min(250, std::max(15, base_target_bits - lp2_bits + 2)));

    for (int poly_i = 0; poly_i < polys_to_run; ++poly_i)
    {
        if ((poly_i & 15) == 0 && stop_requested())
            return found_relations;

        if (poly_i > 0)
        {
            const int gray_curr = poly_i ^ (poly_i >> 1);
            const int gray_prev = (poly_i - 1) ^ ((poly_i - 1) >> 1);
            const int diff = gray_curr ^ gray_prev;
            const int v = static_cast<int>(std::countr_zero(static_cast<uint32_t>(diff))) + 1;

            const int new_sign = -signs[v];
            signs[v] = new_sign;
            const auto &row_v = B_terms_mod_p[v];
            Big two_Bv = B_terms[v];
            two_Bv.shl1();

            if (new_sign == -1)
            {
                b_signed_add(two_Bv, true);
                update_roots(num_primes, cfg.prime_base.data(), row_v.data(), sol1.data(),
                             sol2.data(), true);
            }
            else
            {
                b_signed_add(two_Bv, false);
                update_roots(num_primes, cfg.prime_base.data(), row_v.data(), sol1.data(),
                             sol2.data(), false);
            }
        }

        std::memset(sieve_logs.data(), 0, B_size);
        uint8_t *__restrict s_ptr = sieve_logs.data();

        // 1. Просеивание малых и средних простых.
        for (u32 idx = cfg.small_prime_skip; idx < cfg.split_idx; ++idx)
        {
            if (inv_A_mod_p[idx] == 0)
                continue;
            const u32 p = cfg.prime_base[idx];
            const uint8_t lp = cfg.log_primes[idx];
            const u32 s1 = sol1[idx];
            const u32 s2 = sol2[idx];
            for (u32 s_idx = s1; s_idx < B_size; s_idx += p)
                s_ptr[s_idx] += lp;
            if (s2 != s1)
                for (u32 s_idx = s2; s_idx < B_size; s_idx += p)
                    s_ptr[s_idx] += lp;
        }

        // 2. Просеивание крупных простых (p >= B_size).
        in_block_large1.clear();
        in_block_large2.clear();
        for (u32 idx = cfg.split_idx; idx < num_primes; ++idx)
        {
            const u32 s1 = sol1[idx];
            const u32 s2 = sol2[idx];
            const uint8_t lp = cfg.log_primes[idx];
            if (s1 < B_size)
            {
                s_ptr[s1] += lp;
                in_block_large1.push_back((static_cast<u64>(s1) << 32) | idx);
            }
            if (s2 < B_size)
            {
                s_ptr[s2] += lp;
                in_block_large2.push_back((static_cast<u64>(s2) << 32) | idx);
            }
        }

        // 3. Поиск выживших кандидатов.
        // Скан блока — самая дорогая часть после просева: на полином приходится
        // B_size побайтовых сравнений. Сравниваем по 16 байт (SSE2), а результат
        // сразу превращаем в индексы через movemask. Байты логарифмов знаковые
        // только формально, поэтому перед cmpgt сдвигаем диапазон на 0x80.
        survivor_idx.clear();
        std::fill(survivor_bits.begin(), survivor_bits.end(), 0ULL);
#if defined(__SSE2__)
        {
            const __m128i threshold =
                _mm_set1_epi8(static_cast<char>((target_log - 1) ^ 0x80));
            const __m128i flip = _mm_set1_epi8(static_cast<char>(0x80));
            u32 s_idx = 0;
            for (; s_idx + 16 <= B_size; s_idx += 16)
            {
                __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s_ptr + s_idx));
                v = _mm_xor_si128(v, flip);
                unsigned mask =
                    static_cast<unsigned>(_mm_movemask_epi8(_mm_cmpgt_epi8(v, threshold)));
                while (mask)
                {
                    survivor_idx.push_back(s_idx + static_cast<u32>(std::countr_zero(mask)));
                    mask &= mask - 1;
                }
            }
            for (; s_idx < B_size; ++s_idx)
                if (s_ptr[s_idx] >= target_log)
                    survivor_idx.push_back(s_idx);
        }
#else
        for (u32 s_idx = 0; s_idx < B_size; ++s_idx)
            if (s_ptr[s_idx] >= target_log)
                survivor_idx.push_back(s_idx);
#endif

        if (survivor_idx.empty())
            continue;

        for (u32 s_idx : survivor_idx)
            survivor_bits[s_idx >> 6] |= 1ULL << (s_idx & 63);

        // Крупные простые (p >= B_size) делят Q(x) только там, где их корень
        // попал в блок. Отбираем попадания в выжившие позиции: проверка — один
        // бит из карты в L1 вместо поиска в хеш-таблице на каждое простое.
        // Позиция уже лежит рядом с индексом, поэтому читать sol1/sol2 заново
        // не нужно — проход идёт по плотному массиву.
        hit_pairs.clear();
        for (u64 h : in_block_large1)
        {
            const u32 pos = static_cast<u32>(h >> 32);
            if ((survivor_bits[pos >> 6] >> (pos & 63)) & 1ULL)
                hit_pairs.push_back(h);
        }
        for (u64 h : in_block_large2)
        {
            const u32 pos = static_cast<u32>(h >> 32);
            if ((survivor_bits[pos >> 6] >> (pos & 63)) & 1ULL)
                hit_pairs.push_back(h);
        }

        // 4. Точное разложение выживших кандидатов.
        for (u32 s_idx : survivor_idx)
        {
            large_divs.clear();
            for (u64 h : hit_pairs)
                if (static_cast<u32>(h >> 32) == s_idx)
                    large_divs.push_back(static_cast<u32>(h));

            const u32 mask = half_B ^ s_idx;
            const i64 y = static_cast<i64>(half_B ^ mask) - static_cast<i64>(half_B);

            const bool y_neg = (y < 0);
            const Big ay = A * Big(static_cast<u64>(y_neg ? -y : y));
            Big X_raw;
            bool x_neg;
            if (y_neg == B_neg)
            {
                X_raw = ay + B_mag;
                x_neg = y_neg;
            }
            else if (ay >= B_mag)
            {
                X_raw = ay - B_mag;
                x_neg = y_neg;
            }
            else
            {
                X_raw = B_mag - ay;
                x_neg = B_neg;
            }

            const Big X2 = X_raw * X_raw;
            if (X2 == cfg.kM)
                continue;

            int sign = 1;
            Big diff;
            if (X2 >= cfg.kM)
            {
                diff = X2 - cfg.kM;
                sign = 1;
            }
            else
            {
                diff = cfg.kM - X2;
                sign = -1;
            }

            Big temp = diff / A;
            exps.clear();
            for (u32 a_idx : chosen_indices)
                exps.push_back({a_idx, 1});

            if ((temp.low_word() & 1) == 0)
            {
                u32 cnt2 = 0;
                while (!temp.is_zero() && (temp.low_word() & 1) == 0)
                {
                    const u32 step = (temp.low_word() != 0)
                                         ? std::countr_zero(temp.low_word())
                                         : 64U;
                    temp.shr(step);
                    cnt2 += step;
                }
                if (cnt2 != 0)
                    exps.push_back({0, cnt2});
            }

            for (u32 idx : large_divs)
            {
                const u32 p = cfg.prime_base[idx];
                if (temp.mod_u64(p) == 0)
                {
                    u32 cnt = 0;
                    while (temp.mod_u64(p) == 0)
                    {
                        temp.divrem_u64(p);
                        ++cnt;
                    }
                    merge_exp(exps, idx, cnt);
                }
            }

            for (u32 idx = 1; idx < cfg.split_idx; ++idx)
            {
                const u32 p = cfg.prime_base[idx];
                // Остаток через предвычисленное floor(2^32 / p) вместо деления.
                const u32 q_est = static_cast<u32>((static_cast<u64>(s_idx) * p_recip[idx]) >> 32);
                u32 rem_s = s_idx - q_est * p;
                if (rem_s >= p)
                    rem_s -= p;
                if (rem_s == sol1[idx] || rem_s == sol2[idx] || is_chosen[idx])
                {
                    if (temp.mod_u64(p) == 0)
                    {
                        u32 cnt = 0;
                        while (temp.mod_u64(p) == 0)
                        {
                            temp.divrem_u64(p);
                            ++cnt;
                        }
                        merge_exp(exps, idx, cnt);
                        if (temp.bit_width() <= 1)
                            break;
                    }
                }
            }

            if (temp.bit_width() > 63)
                continue;
            u64 cof = temp.low_word();
            if (cof > cfg.lp2_bound)
                continue;

            if (cof > 1 && cof <= cfg.last_prime)
            {
                auto pit = cfg.prime_to_idx.find(cof);
                if (pit != cfg.prime_to_idx.end())
                {
                    merge_exp(exps, pit->second, 1);
                    cof = 1;
                }
                else
                {
                    continue;
                }
            }

            Big X_val = X_raw % cfg.M;
            if (x_neg && !X_val.is_zero())
                X_val = cfg.M - X_val;

            if (cof == 1)
            {
                found_relations.push_back({1, 1, X_val, sign, exps});
            }
            else if (cof <= cfg.lp_bound)
            {
                found_relations.push_back({1, cof, X_val, sign, exps});
            }
            else if (cof >= cfg.min_2lp_comp)
            {
                if (!is_prime64(cof))
                {
                    u64 f1 = 0, f2 = 0;
                    if (brent_factor_two_primes(cof, f1, f2))
                    {
                        if (f1 > cfg.last_prime && f1 <= cfg.lp_bound
                            && f2 > cfg.last_prime && f2 <= cfg.lp_bound)
                        {
                            if (is_prime64(f1) && is_prime64(f2))
                                found_relations.push_back({f1, f2, X_val, sign, exps});
                        }
                    }
                }
            }
        }
    }
    return found_relations;
}

// --- Автоподбор параметров --------------------------------------------------

struct AutoParams
{
    u32 prime_bound;
    int block_bits;
    int num_a_factors;
};

inline AutoParams auto_tune_siqs(const Big &M)
{
    const int bits = static_cast<int>(M.bit_width());
    const double ln_M = bits * 0.6931471805599453;
    const double ln_ln_M = std::log(ln_M);

    // Граница факторной базы. Свободный член подобран по замерам: на 160/192/205/
    // 224/256 битах оптимум лежит примерно втрое ниже значения, которое давала
    // прежняя формула (+3.3). Меньше база — больше полиномов, но каждый дешевле:
    // при +2.0 время падает в 1.3…2.6 раза, а матрица Гаусса получается меньше
    // (см. таблицу производительности в README).
    const double exponent = 0.357 * std::sqrt(ln_M * ln_ln_M) + 2.0;
    const double raw_bound = std::exp(exponent);

    u32 prime_bound = static_cast<u32>(std::round(raw_bound / 10000.0) * 10000.0);
    const bool base_at_floor = prime_bound < 50000U;
    prime_bound = std::max(50000U, std::min(prime_bound, 8000000U));

    // Ширина блока. На 160–192 битах замеры дают 17 бит (B = 65536) против 16,
    // на 128 — наоборот, 16: слишком широкий интервал раздувает |Q(x)| и
    // соотношений становится меньше, чем экономит просев.
    int block_bits = 16;
    if (bits >= 150)
        block_bits = 17;
    if (bits >= 270)
        block_bits = 18;

    const double half_B = static_cast<double>(1U << (block_bits - 1));
    const double log2_target_A = 0.5 * (bits + 1.0) - std::log2(half_B);
    const double log2_avg_q = std::log2(static_cast<double>(prime_bound)) * 0.62;
    int s = static_cast<int>(std::round(log2_target_A / log2_avg_q));
    // Если база упёрлась в нижний предел, годных q заметно меньше, чем предсказывает
    // формула, и лишний множитель в A только портит полиномы: на 144–176 битах
    // s на единицу меньше даёт 1.2–1.6 раза (например, 160 бит: 0.62 → 0.40 с).
    if (base_at_floor)
        --s;
    // Потолок 10: на 300 битах формула просит 10.7, но 11 множителей дают
    // измеримо худший результат (1546 с против 1508 с) — выигрыш от вдвое более
    // редкой подготовки A не покрывает ухудшения самих полиномов.
    s = std::max(5, std::min(s, 10));

    return {prime_bound, block_bits, s};
}

inline u32 find_best_multiplier(const Big &M, const std::vector<u32> &small_primes)
{
    u32 best_k = 1;
    double best_score = -1e9;
    static const u32 cands[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21,
                                22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 41,
                                42, 43, 46, 47};
    for (u32 k : cands)
    {
        // kM должно помещаться в 512 бит без переполнения.
        if (M.bit_width() + static_cast<u32>(std::bit_width(k)) > 512)
            continue;
        const Big kM = M * Big(k);
        double score = -0.5 * std::log(static_cast<double>(k));
        const u64 mod8 = kM.low_word() & 7ULL;
        if (mod8 == 1)
            score += 2.0 * std::log(2.0);
        else if (mod8 == 5)
            score += 1.0 * std::log(2.0);
        for (size_t i = 1; i < std::min<size_t>(40, small_primes.size()); ++i)
        {
            const u32 p = small_primes[i];
            const u64 rem = kM.mod_u64(p);
            if (rem == 0)
                score += std::log(static_cast<double>(p)) / p;
            else if (pow_mod64(rem, (p - 1) / 2, p) == 1)
                score += 2.0 * std::log(static_cast<double>(p)) / (p - 1);
        }
        if (score > best_score)
        {
            best_score = score;
            best_k = k;
        }
    }
    return best_k;
}

inline std::vector<u32> sieve_primes(u32 limit)
{
    std::vector<uint8_t> is_p(static_cast<size_t>(limit) + 1, 1);
    is_p[0] = 0;
    if (limit >= 1)
        is_p[1] = 0;
    for (u32 i = 2; static_cast<u64>(i) * i <= limit; ++i)
        if (is_p[i])
            for (u64 j = static_cast<u64>(i) * i; j <= limit; j += i)
                is_p[static_cast<size_t>(j)] = 0;
    std::vector<u32> primes;
    for (u32 i = 2; i <= limit; ++i)
        if (is_p[i])
            primes.push_back(i);
    return primes;
}

// --- Основной запуск --------------------------------------------------------

inline std::optional<Big> run_siqs(const Big &M, u32 prime_bound, int block_bits,
                                   int num_a_factors,
                                   const std::function<void(const std::string &)> &status_cb,
                                   const std::function<void(double, const std::string &)> &progress_cb)
{
    const auto t0 = std::chrono::steady_clock::now();
    auto notify = [&](const std::string &msg)
    {
        if (status_cb)
            status_cb(msg);
    };
    auto report = [&](double fraction, const std::string &stage)
    {
        if (progress_cb)
            progress_cb(fraction, stage);
    };

    const std::vector<u32> raw_primes = sieve_primes(prime_bound);
    const u32 k = find_best_multiplier(M, raw_primes);
    const Big kM = M * Big(k);
    report(0.02, "подготовка решета");

    SIQSConfig cfg;
    cfg.M = M;
    cfg.kM = kM;
    cfg.B_size = 1U << block_bits;
    cfg.half_B = cfg.B_size / 2;
    cfg.s = num_a_factors;
    cfg.small_prime_skip = 17;

    for (u32 p : raw_primes)
    {
        if (p == 2)
        {
            cfg.prime_base.push_back(2);
            const u32 r = static_cast<u32>(kM.low_word() & 1ULL);
            cfg.root1.push_back(r);
            cfg.root2.push_back(r);
            cfg.log_primes.push_back(1);
        }
        else
        {
            const u64 rem = kM.mod_u64(p);
            if (rem == 0 || pow_mod64(rem, (p - 1) / 2, p) == 1)
            {
                const u32 r = tonelli_shanks(rem, p);
                cfg.prime_base.push_back(p);
                cfg.root1.push_back(r);
                cfg.root2.push_back((r != p - r) ? (p - r) : r);
                cfg.log_primes.push_back(
                    static_cast<uint8_t>(std::llround(std::log2(static_cast<double>(p)))));
            }
        }
    }

    cfg.num_primes = static_cast<u32>(cfg.prime_base.size());
    if (cfg.num_primes < 64)
        return std::nullopt;

    for (u32 i = 0; i < cfg.num_primes; ++i)
        cfg.prime_to_idx[cfg.prime_base[i]] = i;
    cfg.last_prime = cfg.prime_base.back();
    cfg.min_2lp_comp = cfg.last_prime * cfg.last_prime;
    cfg.lp_bound = static_cast<u64>(prime_bound) * 30ULL;
    cfg.lp2_bound = cfg.lp_bound * cfg.lp_bound;

    cfg.split_idx = cfg.num_primes;
    for (u32 i = 0; i < cfg.num_primes; ++i)
    {
        if (cfg.prime_base[i] >= cfg.B_size)
        {
            cfg.split_idx = i;
            break;
        }
    }

    UInt<1024> two_kM(kM);
    two_kM.shl(1);
    UInt<1024> target_A = ops::isqrt<1024>(two_kM);
    target_A.divrem_u64(cfg.half_B);
    const double log2_A = static_cast<double>(target_A.bit_width());
    const u32 ideal_q = static_cast<u32>(
        std::max(50.0, std::pow(2.0, log2_A / static_cast<double>(cfg.s))));

    u32 center_idx = cfg.num_primes / 2;
    for (u32 i = 0; i < cfg.num_primes; ++i)
    {
        if (cfg.prime_base[i] >= ideal_q)
        {
            center_idx = i;
            break;
        }
    }
    u32 pool_start = 15;
    if (cfg.num_primes > 200)
    {
        const u32 lo = (center_idx > 60) ? (center_idx - 60) : 15U;
        pool_start = std::max(15U, std::min(lo, cfg.num_primes - 130));
    }
    for (u32 i = pool_start; i < std::min(cfg.num_primes - 1, pool_start + 120); ++i)
        cfg.pool_indices.push_back(i);

    unsigned threads = cpu::worker_thread_count();
    if (threads == 0)
        threads = 1;

    notify("SIQS: " + std::to_string(M.bit_width()) + " бит, множитель k="
         + std::to_string(k) + ", база: " + std::to_string(cfg.num_primes)
         + " простых до " + std::to_string(prime_bound)
         + ", полиномов на A: " + std::to_string(1ULL << (cfg.s - 1))
         + ", потоков: " + std::to_string(threads));
    report(0.05, "сбор соотношений");

    MarkowitzGaussianGF2 gauss;
    gauss.status_cb = status_cb;
    gauss.progress_cb = [&](double f, const std::string &t)
    {
        report(0.86 + 0.12 * f, t);
    };
    LargePrimeGraph lp_graph(M);
    u64 total_polys = 0;
    u64 batch_seed = 1000;
    const int polys_per_A = static_cast<int>(1ULL << (cfg.s - 1));
    // Ориентир для полосы прогресса: для готовности матрицы нужно примерно
    // столько соотношений, сколько активных столбцов (простых базы).
    const double relation_target = static_cast<double>(std::max<u32>(cfg.num_primes, 64));

    // Нарастающая порция полиномов на поток: начинаем с двух и удваиваем, пока
    // матрица не готова. Так мелкие числа не собирают соотношения «про запас»
    // (иначе Гаусс считается на матрице в десятки раз больше нужной), а крупные
    // быстро выходят на полный батч.
    int poly_limit = std::min(polys_per_A, 2);

    while (true)
    {
        if (stop_requested())
            return std::nullopt;

        std::vector<std::vector<RawRelation>> batch_results(threads);

        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (unsigned t = 0; t < threads; ++t)
        {
            pool.emplace_back([&, t]()
            {
                cpu::pin_current_thread(t);
                batch_results[t] = process_one_A_batch(cfg, batch_seed + t + 1, poly_limit);
            });
        }
        for (auto &th : pool)
            if (th.joinable())
                th.join();
        batch_seed += threads;
        total_polys += static_cast<u64>(threads) * static_cast<u64>(poly_limit);
        poly_limit = std::min(polys_per_A, poly_limit * 2);

        for (unsigned t = 0; t < threads; ++t)
        {
            for (const auto &r : batch_results[t])
            {
                FullRelation full_rel;
                if (lp_graph.add_edge_and_extract_cycle(r, full_rel))
                    gauss.add_relation(std::move(full_rel));
            }
        }

        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - t0).count();
        notify("SIQS: " + std::to_string(static_cast<long long>(elapsed)) + " с, полиномов "
             + std::to_string(total_polys) + ", соотношений: "
             + std::to_string(gauss.relations.size()) + ", простых (>=2 вх.): "
             + std::to_string(gauss.multi_hit_primes) + " (синглетонов: "
             + std::to_string(gauss.singleton_primes) + ", циклов 2LP: "
             + std::to_string(lp_graph.cycles_found) + ")");

        // Полоса сбора соотношений: 5%..85% по доле от ориентира. Если матрица
        // ещё не готова, полоса просто упирается в верхнюю границу — дальнейшие
        // батчи продолжают уточнять матрицу, а не «съедают» прогресс.
        report(0.05 + 0.80 * std::min(1.0, static_cast<double>(gauss.relations.size()) / relation_target),
               "сбор соотношений");

        if (!gauss.ready_to_check_rank())
            continue;

        report(0.86, "решение матрицы GF(2)");
        const auto tg0 = std::chrono::steady_clock::now();
        const auto deps = gauss.solve_pruned_matrix();
        for (const auto &dep : deps)
        {
            if (stop_requested())
                return std::nullopt;

            Big LHS(1);
            std::unordered_map<u32, u32> total_exps;
            for (u32 idx_r : dep)
            {
                LHS = Big::mul_mod(LHS, gauss.relations[idx_r].X_val, M);
                for (const auto &fe : gauss.relations[idx_r].exps)
                    total_exps[fe.idx] += fe.cnt;
            }

            Big RHS(1);
            for (const auto &kv : total_exps)
            {
                const u32 p_val = cfg.prime_base[kv.first];
                const Big p_pow = pow_mod_big(Big(p_val), kv.second / 2, M);
                RHS = Big::mul_mod(RHS, p_pow, M);
            }

            const Big diff = (LHS >= RHS) ? (LHS - RHS) : (RHS - LHS);
            Big sum = LHS + RHS;
            if (sum >= M)
                sum.sub(M);

            const Big candidates[2] = {diff, sum};
            for (const Big &cand : candidates)
            {
                const Big g = ops::gcd(cand, M);
                if (g > Big(1) && g < M)
                {
                    const auto t1 = std::chrono::steady_clock::now();
                    const double total_s = std::chrono::duration<double>(t1 - t0).count();
                    const double gauss_s = std::chrono::duration<double>(t1 - tg0).count();
                    notify("SIQS: успех за " + std::to_string(total_s) + " с (Гаусс: "
                         + std::to_string(gauss_s) + " с), полиномов: "
                         + std::to_string(total_polys) + ", делитель: " + g.toString());
                    report(1.0, "делитель найден");
                    return g;
                }
            }
        }
    }
}

// Точка входа для конвейера: факторизует n, возвращает нетривиальный делитель.
// progress_cb (необязательный) получает долю 0..1 внутри решета и название
// текущей фазы — конвейер отображает её в общую полосу прогресса.
template <typename T>
std::optional<T> factor(const T &n,
                        std::function<void(const std::string &)> status_cb = nullptr,
                        std::function<void(double, const std::string &)> progress_cb = nullptr)
{
    Big M(n);
    if (M <= Big(3))
        return std::nullopt;
    if ((M.low_word() & 1ULL) == 0)
        return T(2);

    const AutoParams params = auto_tune_siqs(M);
    if (status_cb)
        status_cb("SIQS: автоподбор: prime_bound=" + std::to_string(params.prime_bound)
                  + ", block_bits=" + std::to_string(params.block_bits)
                  + ", num_a_factors=" + std::to_string(params.num_a_factors));

    const auto divisor = run_siqs(M, params.prime_bound, params.block_bits,
                                  params.num_a_factors, status_cb, progress_cb);
    if (!divisor.has_value())
        return std::nullopt;
    return T(*divisor);
}

} // namespace bignum::siqs_xor
