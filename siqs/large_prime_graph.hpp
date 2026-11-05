#pragma once

// ---------------------------------------------------------------------------
// siqs/large_prime_graph.hpp — граф циклов 2LP.
//
// Ориентированный остовный лес по большим простым: каждое соотношение с
// двумя большими простыми — ребро, цикл в лесу даёт готовую строку матрицы.
// ---------------------------------------------------------------------------

#include "common.hpp"
#include "relation.hpp"

namespace bignum::siqs_xor {

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
        // Совместимость с 2LP-логикой ниже: полное соотношение — (1,1),
        // одно большое простое — (1,p), два — (p,q). 3LP пока не сюда.
        u64 u = 1, v = 1;
        if (rel.n_lp >= 1)
            v = rel.lp[0];
        if (rel.n_lp >= 2)
        {
            u = rel.lp[0];
            v = rel.lp[1];
        }
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

} // namespace bignum::siqs_xor
