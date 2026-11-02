// Компиляция (GCC / Clang):
// g++ -O3 -march=native -fopenmp siqs_xor.cpp -o siqs_xor
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <chrono>
#include <omp.h>

using u32 = uint32_t;
using u64 = uint64_t;
using i64 = int64_t;
using u128 = unsigned __int128;

// Фиксированное 512-битное беззнаковое целое на стеке (8 x 64 бит)
struct U512
{
    u64 w[8];

    U512() { std::memset(w, 0, sizeof(w)); }
    U512(u64 v)
    {
        std::memset(w, 0, sizeof(w));
        w[0] = v;
    }

    bool is_zero() const
    {
        for (int i = 0; i < 8; ++i)
            if (w[i] != 0)
                return false;
        return true;
    }

    int cmp(const U512 &o) const
    {
        for (int i = 7; i >= 0; --i)
        {
            if (w[i] < o.w[i])
                return -1;
            if (w[i] > o.w[i])
                return 1;
        }
        return 0;
    }
    bool operator<(const U512 &o) const { return cmp(o) < 0; }
    bool operator<=(const U512 &o) const { return cmp(o) <= 0; }
    bool operator>(const U512 &o) const { return cmp(o) > 0; }
    bool operator>=(const U512 &o) const { return cmp(o) >= 0; }
    bool operator==(const U512 &o) const { return cmp(o) == 0; }
    bool operator!=(const U512 &o) const { return cmp(o) != 0; }

    int bit_length() const
    {
        for (int i = 7; i >= 0; --i)
        {
            if (w[i] != 0)
                return i * 64 + (64 - __builtin_clzll(w[i]));
        }
        return 0;
    }

    bool get_bit(int b) const
    {
        return (w[b >> 6] >> (b & 63)) & 1ULL;
    }
    void set_bit(int b)
    {
        w[b >> 6] |= (1ULL << (b & 63));
    }

    void add(const U512 &o)
    {
        u64 carry = 0;
        for (int i = 0; i < 8; ++i)
        {
            u128 s = (u128)w[i] + o.w[i] + carry;
            w[i] = (u64)s;
            carry = (u64)(s >> 64);
        }
    }

    void sub(const U512 &o)
    {
        u64 borrow = 0;
        for (int i = 0; i < 8; ++i)
        {
            u128 d = (u128)w[i] - o.w[i] - borrow;
            w[i] = (u64)d;
            borrow = (d >> 127) ? 1 : 0;
        }
    }

    void shl1()
    {
        u64 c = 0;
        for (int i = 0; i < 8; ++i)
        {
            u64 next_c = w[i] >> 63;
            w[i] = (w[i] << 1) | c;
            c = next_c;
        }
    }

    void shr_k(int k)
    {
        while (k > 0)
        {
            int step = std::min(k, 63);
            u64 c = 0;
            for (int i = 7; i >= 0; --i)
            {
                u64 next_c = w[i] << (64 - step);
                w[i] = (w[i] >> step) | c;
                c = next_c;
            }
            k -= step;
        }
    }

    // Сверхбыстрое деление с остатком на 64-битное число (для факторной базы!)
    u64 divrem_u64(u64 d)
    {
        u128 rem = 0;
        for (int i = 7; i >= 0; --i)
        {
            if (w[i] == 0 && rem == 0)
                continue;
            u128 cur = (rem << 64) | w[i];
            w[i] = (u64)(cur / d);
            rem = cur % d;
        }
        return (u64)rem;
    }

    u64 mod_u64(u64 d) const
    {
        u128 rem = 0;
        for (int i = 7; i >= 0; --i)
        {
            if (w[i] == 0 && rem == 0)
                continue;
            rem = ((rem << 64) | w[i]) % d;
        }
        return (u64)rem;
    }
};

// Вспомогательные функции для U512
U512 add512(U512 a, const U512 &b)
{
    a.add(b);
    return a;
}
U512 sub512(U512 a, const U512 &b)
{
    a.sub(b);
    return a;
}

U512 mul512(const U512 &a, const U512 &b)
{
    U512 res;
    for (int i = 0; i < 8; ++i)
    {
        if (a.w[i] == 0)
            continue;
        u64 carry = 0;
        for (int j = 0; i + j < 8; ++j)
        {
            u128 cur = (u128)a.w[i] * b.w[j] + res.w[i + j] + carry;
            res.w[i + j] = (u64)cur;
            carry = (u64)(cur >> 64);
        }
    }
    return res;
}

U512 mul_u64(const U512 &a, u64 b)
{
    U512 res;
    u64 carry = 0;
    for (int i = 0; i < 8; ++i)
    {
        u128 cur = (u128)a.w[i] * b + carry;
        res.w[i] = (u64)cur;
        carry = (u64)(cur >> 64);
    }
    return res;
}

void divmod512(const U512 &num, const U512 &den, U512 &q, U512 &r)
{
    q = U512(0);
    r = U512(0);
    int bits = num.bit_length();
    for (int i = bits - 1; i >= 0; --i)
    {
        r.shl1();
        if (num.get_bit(i))
            r.w[0] |= 1ULL;
        if (r >= den)
        {
            r.sub(den);
            q.set_bit(i);
        }
    }
}

U512 div512(const U512 &a, const U512 &b)
{
    U512 q, r;
    divmod512(a, b, q, r);
    return q;
}
U512 mod512(const U512 &a, const U512 &b)
{
    U512 q, r;
    divmod512(a, b, q, r);
    return r;
}

// Умножение по модулю M без переполнения (через бинарное удвоение, 100% точность до 500 бит)
U512 mul_mod512(U512 a, U512 b, const U512 &m)
{
    a = mod512(a, m);
    b = mod512(b, m);
    U512 res(0);
    int bits = b.bit_length();
    for (int i = 0; i < bits; ++i)
    {
        if (b.get_bit(i))
        {
            res.add(a);
            if (res >= m)
                res.sub(m);
        }
        a.shl1();
        if (a >= m)
            a.sub(m);
    }
    return res;
}

U512 pow_mod512(U512 base, u64 exp, const U512 &m)
{
    U512 res(1);
    base = mod512(base, m);
    while (exp > 0)
    {
        if (exp & 1)
            res = mul_mod512(res, base, m);
        base = mul_mod512(base, base, m);
        exp >>= 1;
    }
    return res;
}

U512 gcd512(U512 a, U512 b)
{
    while (!b.is_zero())
    {
        U512 r = mod512(a, b);
        a = b;
        b = r;
    }
    return a;
}

// Обратный элемент a^(-1) mod m (100% точный модулярный Евклид без знаковых флагов!)
U512 inv_mod512(U512 a, const U512& m) {
    a = mod512(a, m);
    U512 u = a, v = m;
    U512 x1(1), x2(0);
    while (!u.is_zero()) {
        U512 q, r;
        divmod512(v, u, q, r);
        v = u;
        u = r;
        U512 qx1 = mul_mod512(q, x1, m);
        // Вычисляем nx = (x2 - qx1) mod m чисто в положительных вычетах [0 .. m-1]:
        U512 nx = (x2 >= qx1) ? sub512(x2, qx1) : sub512(m, sub512(qx1, x2));
        x2 = x1;
        x1 = nx;
    }
    return x2;
}

U512 isqrt512(const U512 &n)
{
    if (n.is_zero())
        return U512(0);
    U512 res(0);
    int b = ((n.bit_length() - 1) >> 1);
    for (int i = b; i >= 0; --i)
    {
        U512 cand = res;
        cand.set_bit(i);
        if (mul512(cand, cand) <= n)
            res = cand;
    }
    return res;
}

U512 from_dec_string(const std::string &s)
{
    U512 res(0);
    for (char c : s)
    {
        if (c >= '0' && c <= '9')
        {
            res = mul_u64(res, 10);
            res.add(U512((u64)(c - '0')));
        }
    }
    return res;
}

std::string to_dec_string(U512 v)
{
    if (v.is_zero())
        return "0";
    std::string s;
    while (!v.is_zero())
    {
        u64 rem = v.divrem_u64(10);
        s.push_back((char)('0' + rem));
    }
    std::reverse(s.begin(), s.end());
    return s;
}

// Быстрая 64-битная модулярная арифметика через аппаратный __int128
inline u64 mul_mod64(u64 a, u64 b, u64 m)
{
    return (u64)(((u128)a * b) % m);
}

inline u64 pow_mod64(u64 base, u64 exp, u64 m)
{
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

u64 inv_mod64(u64 a, u64 m)
{
    i64 t = 0, newt = 1;
    i64 r = (i64)m, newr = (i64)(a % m);
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
        t += (i64)m;
    return (u64)t;
}

// Детерминированный тест Миллера — Рабина для u64
bool is_prime64(u64 n)
{
    if (n < 2)
        return false;
    if (n == 2 || n == 3 || n == 5 || n == 7)
        return true;
    if ((n & 1) == 0 || n % 3 == 0)
        return false;
    u64 d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;
    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17};
    for (u64 a : bases)
    {
        if (n <= a)
            break;
        u64 x = pow_mod64(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool comp = true;
        for (int i = 0; i < r - 1; ++i)
        {
            x = mul_mod64(x, x, n);
            if (x == n - 1)
            {
                comp = false;
                break;
            }
        }
        if (comp)
            return false;
    }
    return true;
}

// Сверхбыстрый метод Брента — Полларда на регистрах CPU
bool brent_factor_two_primes(u64 n, u64 &f1, u64 &f2)
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
                u64 lim = std::min(m, r - k);
                for (u64 i = 0; i < lim; ++i)
                {
                    y = (mul_mod64(y, y, n) + c) % n;
                    u64 diff = (x > y) ? (x - y) : (y - x);
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
                u64 diff = (x > ys) ? (x - ys) : (ys - x);
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

// Алгоритм Тонелли — Шэнкса: r^2 ≡ n (mod p)
u32 tonelli_shanks(u64 n, u32 p)
{
    n %= p;
    if (n == 0)
        return 0;
    if (p == 2)
        return (u32)n;
    if ((p & 3) == 3)
        return (u32)pow_mod64(n, (p + 1) / 4, p);
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
        u64 b = pow_mod64(c, 1ULL << (m - i - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }
    return (u32)r;
}

constexpr int Z0 = 0;

struct FactorExp
{
    u32 idx;
    u32 cnt;
};

struct RawRelation
{
    u64 u_lp;
    u64 v_lp;
    U512 X_val;
    int sign;
    std::vector<FactorExp> exps;
};

struct FullRelation
{
    U512 X_val;
    int sign;
    std::vector<FactorExp> exps;
};

static void merge_exps(std::vector<FactorExp> &dst, const std::vector<FactorExp> &src)
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

// Сверхбыстрый граф циклов 2LP на ориентированном остовном лесе
class LargePrimeGraph
{
    struct TreeEdge
    {
        u64 parent_node;
        U512 X_val;
        int sign;
        std::vector<FactorExp> exps;
        u64 p1, p2;
    };

    U512 M;
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
        u64 root = curr;
        curr = v;
        while (true)
        {
            auto it = uf_parent.find(curr);
            if (it == uf_parent.end() || it->second == root)
                break;
            u64 nxt = it->second;
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
        for (int i = (int)path.size() - 1; i >= 0; --i)
        {
            u64 child = path[i].first;
            TreeEdge te = path[i].second;
            u64 par = te.parent_node;
            te.parent_node = child;
            tree_parent[par] = std::move(te);
        }
        tree_parent.erase(node);
    }

public:
    u64 num_edges = 0;
    u64 cycles_found = 0;

    LargePrimeGraph(const U512 &m_val) : M(m_val) {}

    bool add_edge_and_extract_cycle(const RawRelation &rel, FullRelation &out_rel)
    {
        u64 u = rel.u_lp, v = rel.v_lp;
        if (u == 1 && v == 1)
        {
            out_rel = {rel.X_val, rel.sign, rel.exps};
            return true;
        }
        if (u == v)
        {
            U512 inv_u = inv_mod512(U512(u), M);
            out_rel = {mul_mod512(rel.X_val, inv_u, M), rel.sign, rel.exps};
            return true;
        }

        u64 ru = find_root(u);
        u64 rv = find_root(v);
        u64 ep1 = (u > 1) ? u : 0;
        u64 ep2 = (v > 1) ? v : 0;

        if (ru != rv)
        {
            reroot_tree(u);
            tree_parent[u] = {v, rel.X_val, rel.sign, rel.exps, ep1, ep2};
            uf_parent[ru] = rv;
            num_edges++;
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
        {
            curr = tree_parent[curr].parent_node;
        }
        u64 lca = curr;

        U512 comb_X = rel.X_val;
        int comb_sign = rel.sign;
        std::vector<FactorExp> comb_exps = rel.exps;
        std::unordered_map<u64, u32> lp_counts;
        if (ep1)
            lp_counts[ep1]++;
        if (ep2)
            lp_counts[ep2]++;

        for (u64 start_node : {u, v})
        {
            curr = start_node;
            while (curr != lca)
            {
                const auto &te = tree_parent[curr];
                comb_X = mul_mod512(comb_X, te.X_val, M);
                comb_sign *= te.sign;
                merge_exps(comb_exps, te.exps);
                if (te.p1)
                    lp_counts[te.p1]++;
                if (te.p2)
                    lp_counts[te.p2]++;
                curr = te.parent_node;
            }
        }

        U512 denom(1);
        for (const auto &kv : lp_counts)
        {
            U512 lp_pow = pow_mod512(U512(kv.first), kv.second / 2, M);
            denom = mul_mod512(denom, lp_pow, M);
        }
        U512 inv_denom = inv_mod512(denom, M);
        comb_X = mul_mod512(comb_X, inv_denom, M);

        cycles_found++;
        out_rel = {comb_X, comb_sign, std::move(comb_exps)};
        return true;
    }
};

// Решатель Гаусса Марковица на 64-битных машинных словах (SGE + Bitpacked GF(2))
class MarkowitzGaussianGF2
{
public:
    std::vector<FullRelation> relations;
    std::unordered_map<u32, u32> col_freq;
    u32 multi_hit_primes = 0;
    u32 singleton_primes = 0;
    size_t last_check_count = 0;

    void add_relation(FullRelation &&rel)
    {
        for (const auto &fe : rel.exps)
        {
            if (fe.cnt & 1)
            {
                u32 prev = col_freq[fe.idx]++;
                if (prev == 0)
                    singleton_primes++;
                else if (prev == 1)
                {
                    singleton_primes--;
                    multi_hit_primes++;
                }
            }
        }
        relations.push_back(std::move(rel));
    }

    bool ready_to_check_rank() const
    {
        size_t n_rels = relations.size();
        if (n_rels - last_check_count < 500)
            return false;
        size_t total_active = (size_t)multi_hit_primes + singleton_primes;
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

        for (u32 r_idx = 0; r_idx < (u32)relations.size(); ++r_idx)
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
            int s_bit = (relations[r_idx].sign == -1) ? 1 : 0;
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
            {
                if (kv.second.size() == 1)
                    q.push_back(kv.first);
            }
            while (!q.empty())
            {
                u32 c = q.back();
                q.pop_back();
                auto cit = col_to_rows.find(c);
                if (cit != col_to_rows.end() && cit->second.size() == 1)
                {
                    u32 r_del = *cit->second.begin();
                    std::vector<u32> aff(active_rows[r_del].cols.begin(), active_rows[r_del].cols.end());
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
        size_t orig_core = active_rows.size();

        // Сжатие SGE для весов 2..5
        for (size_t max_w : {2, 2, 3, 3, 4, 5})
        {
            std::vector<u32> cands;
            for (const auto &kv : col_to_rows)
            {
                if (kv.second.size() >= 2 && kv.second.size() <= max_w)
                    cands.push_back(kv.first);
            }
            for (u32 c : cands)
            {
                auto cit = col_to_rows.find(c);
                if (cit == col_to_rows.end() || cit->second.size() < 2 || cit->second.size() > max_w)
                    continue;
                if (active_rows.size() <= col_to_rows.size() + 24)
                    break;
                std::vector<u32> r_list(cit->second.begin(), cit->second.end());
                std::sort(r_list.begin(), r_list.end(), [&](u32 a, u32 b)
                          { return active_rows[a].cols.size() < active_rows[b].cols.size(); });
                u32 pivot_r = r_list[Z0];
                if (active_rows[pivot_r].cols.size() > 40)
                    continue;
                for (size_t i = 1; i < r_list.size(); ++i)
                    xor_row_into(pivot_r, r_list[i]);
                remove_row(pivot_r);
            }
            prune_singletons();
        }
        std::cout << "    [Гаусс-триггер C++] Ядро после синглетонов: " << orig_core
                  << " -> сжато SGE до " << active_rows.size()
                  << " строк и " << col_to_rows.size() << " столбцов!\n";
        std::vector<u32> sorted_cols;
        sorted_cols.reserve(col_to_rows.size());
        for (const auto &kv : col_to_rows)
            sorted_cols.push_back(kv.first);
        std::sort(sorted_cols.begin(), sorted_cols.end(), [&](u32 a, u32 b)
                  { return col_to_rows[a].size() > col_to_rows[b].size(); });
        std::unordered_map<u32, u32> col_map;
        for (u32 i = 0; i < (u32)sorted_cols.size(); ++i)
            col_map[sorted_cols[i]] = i;
        u32 sign_bit = (u32)sorted_cols.size();
        u32 num_bits = sign_bit + 1;
        size_t vec_words = (num_bits + 63) / 64;
        std::vector<u32> row_order;
        for (const auto &kv : active_rows)
            row_order.push_back(kv.first);
        std::sort(row_order.begin(), row_order.end(), [&](u32 a, u32 b)
                  { return active_rows[a].cols.size() < active_rows[b].cols.size(); });
        // Упаковываем историю компактно (только по выжившим строкам!)
        size_t hist_words = (row_order.size() + 63) / 64;
        struct PivotEntry
        {
            std::vector<u64> vec;
            std::vector<u64> hist;
        };
        std::unordered_map<int, PivotEntry> pivot_map;
        std::vector<std::vector<u32>> deps;
        for (size_t ord_i = 0; ord_i < row_order.size(); ++ord_i)
        {
            u32 r_idx = row_order[ord_i];
            const auto &ar = active_rows[r_idx];
            std::vector<u64> vec(vec_words, 0);
            for (u32 c : ar.cols)
            {
                auto mit = col_map.find(c);
                if (mit != col_map.end())
                {
                    u32 bit = mit->second;
                    vec[bit >> 6] |= (1ULL << (bit & 63));
                }
            }
            if (ar.s_bit)
                vec[sign_bit >> 6] |= (1ULL << (sign_bit & 63));
            std::vector<u64> hist(hist_words, 0);
            hist[ord_i >> 6] |= (1ULL << (ord_i & 63));
            int msb = -1;
            for (int w = (int)vec_words - 1; w >= 0; --w)
            {
                if (vec[w] != 0)
                {
                    msb = w * 64 + (63 - __builtin_clzll(vec[w]));
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
                    int top_w = msb >> 6;
                    for (int w = 0; w <= top_w; ++w)
                        vec[w] ^= pv[w];
                    for (size_t w = 0; w < hist_words; ++w)
                        hist[w] ^= ph[w];
                    msb = -1;
                    for (int w = top_w; w >= 0; --w)
                    {
                        if (vec[w] != 0)
                        {
                            msb = w * 64 + (63 - __builtin_clzll(vec[w]));
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
                // Восстанавливаем исходные индексы строк через симметрическую разность множеств
                std::vector<u32> comb;
                for (size_t w = 0; w < hist_words; ++w)
                {
                    u64 word = hist[w];
                    while (word)
                    {
                        int b = __builtin_ctzll(word);
                        size_t idx_in_ord = w * 64 + b;
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

// Блок 5 из 6: Многопоточный воркер SIQS (L1-кэш сито + Hit-Map + XOR-маски)
// Весь массив сита `sieve_logs` — это плоский `std::vector<uint8_t>` в кэше процессора, а проверка пустых блоков идет по 8 байт за такт через `uint64_t*`!

struct SIQSConfig
{
    U512 M;
    U512 kM;
    std::vector<u32> prime_base;
    std::unordered_map<u64, u32> prime_to_idx;
    std::vector<u32> root1;
    std::vector<u32> root2;
    std::vector<uint8_t> log_primes;
    u32 num_primes;
    u64 last_prime;
    u64 min_2lp_comp;
    u64 lp_bound;
    u64 lp2_bound;
    u32 B_size;
    u32 half_B;
    int s;
    std::vector<u32> pool_indices;
    u32 small_prime_skip;
    u32 split_idx;
};

std::vector<RawRelation> process_one_A_batch_cpp(const SIQSConfig &cfg, u64 seed)
{
    u32 num_primes = cfg.num_primes;
    u32 B_size = cfg.B_size;
    u32 half_B = cfg.half_B;
    int s = cfg.s;

    std::mt19937_64 rng(seed);
    std::vector<u32> chosen_indices;
    std::vector<u32> pool_copy = cfg.pool_indices;
    for (int i = 0; i < s; ++i)
    {
        std::uniform_int_distribution<size_t> dist(i, pool_copy.size() - 1);
        size_t pick = dist(rng);
        std::swap(pool_copy[i], pool_copy[pick]);
        chosen_indices.push_back(pool_copy[i]);
    }
    std::sort(chosen_indices.begin(), chosen_indices.end());

    std::vector<uint8_t> is_chosen(num_primes, 0);
    U512 A(1);
    for (u32 idx : chosen_indices)
    {
        is_chosen[idx] = 1;
        A = mul_u64(A, cfg.prime_base[idx]);
    }

    std::vector<U512> B_terms(s);
    for (int j = 0; j < s; ++j)
    {
        u32 idx = chosen_indices[j];
        u64 q_j = cfg.prime_base[idx];
        u64 r_j = cfg.root1[idx];
        U512 A_div_q = A;
        A_div_q.divrem_u64(q_j);
        u64 rem_a = A_div_q.mod_u64(q_j);
        u64 inv = inv_mod64(rem_a, q_j);
        u64 gamma = mul_mod64(r_j, inv, q_j);
        if (gamma > q_j / 2)
            gamma = q_j - gamma;
        B_terms[j] = mul_u64(A_div_q, gamma);
    }

    // B = sum(+-B_terms[j]) — знаковая величина (как обычный int в Python-эталоне).
    // Держим модуль и знак: иначе при отрицательном B вычитание в U512 уходит в 2^512.
    U512 B_mag(0);
    bool B_neg = false;
    for (int j = 0; j < s; ++j)
        B_mag.add(B_terms[j]);
    auto b_signed_add = [&](const U512 &term, bool term_neg)
    {
        if (B_neg == term_neg)
        {
            B_mag = add512(B_mag, term);
        }
        else if (B_mag >= term)
        {
            B_mag = sub512(B_mag, term);
        }
        else
        {
            B_mag = sub512(term, B_mag);
            B_neg = term_neg;
        }
        if (B_mag.is_zero())
            B_neg = false;
    };

    std::vector<u32> inv_A_mod_p(num_primes, 0);
    std::vector<std::vector<u32>> B_terms_mod_p(s, std::vector<u32>(num_primes, 0));
    std::vector<u32> sol1(num_primes, 0);
    std::vector<u32> sol2(num_primes, 0);

    for (u32 idx = 1; idx < num_primes; ++idx)
    {
        if (is_chosen[idx])
            continue;
        u32 p = cfg.prime_base[idx];
        u64 a_mod = A.mod_u64(p);
        u32 inv_A = (u32)inv_mod64(a_mod, p);
        inv_A_mod_p[idx] = inv_A;
        for (int j = 0; j < s; ++j)
        {
            u64 b_mod = B_terms[j].mod_u64(p);
            B_terms_mod_p[j][idx] = (u32)((2ULL * b_mod * inv_A) % p);
        }

        u64 B_mod = B_mag.mod_u64(p);
        if (B_neg && B_mod != 0)
            B_mod = (u64)p - B_mod;
        u64 r1 = cfg.root1[idx];
        u64 r2 = cfg.root2[idx];
        // ВАЖНО: считаем в 64-битном u64, чтобы (p + r1 - B_mod) * inv_A не переполняло u32!
        sol1[idx] = (u32)((((p + r1 - B_mod) % p) * (u64)inv_A + (u64)half_B) % p);
        sol2[idx] = (u32)((((p + r2 - B_mod) % p) * (u64)inv_A + (u64)half_B) % p);
    }

    std::vector<int> signs(s, 1);
    int polys_per_A = 1 << (s - 1);
    std::vector<RawRelation> found_relations;

    std::vector<uint8_t> sieve_logs(B_size);
    std::vector<u32> in_block_large1;
    std::vector<u32> in_block_large2;
    in_block_large1.reserve(4096);
    in_block_large2.reserve(4096);

    int base_target_bits = (int)(cfg.kM.bit_length() / 2 + __builtin_ctz(half_B));
    int lp2_bits = 64 - __builtin_clzll(cfg.lp2_bound);
    uint8_t target_log = (uint8_t)std::min(250, std::max(15, base_target_bits - lp2_bits + 2));

    for (int poly_i = 0; poly_i < polys_per_A; ++poly_i)
    {
        if (poly_i > 0) {
            int gray_curr = poly_i ^ (poly_i >> 1);
            int gray_prev = (poly_i - 1) ^ ((poly_i - 1) >> 1);
            int diff = gray_curr ^ gray_prev;
            int v = __builtin_ctz(diff) + 1; // Точный номер изменившегося бита (от 1 до s-1)!

            int new_sign = -signs[v];
            signs[v] = new_sign;
            const auto &row_v = B_terms_mod_p[v];
            U512 two_Bv = B_terms[v];
            two_Bv.shl1();

            if (new_sign == -1)
            {
                b_signed_add(two_Bv, true);
                for (u32 idx = 1; idx < num_primes; ++idx)
                {
                    if (inv_A_mod_p[idx])
                    {
                        u32 p = cfg.prime_base[idx];
                        u32 delta = row_v[idx];
                        u32 s1_n = sol1[idx] + delta;
                        sol1[idx] = (s1_n >= p) ? (s1_n - p) : s1_n;
                        u32 s2_n = sol2[idx] + delta;
                        sol2[idx] = (s2_n >= p) ? (s2_n - p) : s2_n;
                    }
                }
            }
            else
            {
                b_signed_add(two_Bv, false);
                for (u32 idx = 1; idx < num_primes; ++idx)
                {
                    if (inv_A_mod_p[idx])
                    {
                        u32 p = cfg.prime_base[idx];
                        u32 delta = row_v[idx];
                        sol1[idx] = (sol1[idx] >= delta) ? (sol1[idx] - delta) : (sol1[idx] + p - delta);
                        sol2[idx] = (sol2[idx] >= delta) ? (sol2[idx] - delta) : (sol2[idx] + p - delta);
                    }
                }
            }
        }

        std::memset(sieve_logs.data(), 0, B_size);
        uint8_t *__restrict s_ptr = sieve_logs.data();

        // 1. Сверхбыстрое просеивание малых и средних простых в L1/L2 кэше
        for (u32 idx = cfg.small_prime_skip; idx < cfg.split_idx; ++idx)
        {
            if (inv_A_mod_p[idx] == 0)
                continue;
            u32 p = cfg.prime_base[idx];
            uint8_t lp = cfg.log_primes[idx];
            u32 s1 = sol1[idx];
            u32 s2 = sol2[idx];
            for (u32 s_idx = s1; s_idx < B_size; s_idx += p)
                s_ptr[s_idx] += lp;
            if (s2 != s1)
            {
                for (u32 s_idx = s2; s_idx < B_size; s_idx += p)
                    s_ptr[s_idx] += lp;
            }
        }

        // 2. Просеивание крупных простых (p >= B_size)
        in_block_large1.clear();
        in_block_large2.clear();
        for (u32 idx = cfg.split_idx; idx < num_primes; ++idx)
        {
            u32 s1 = sol1[idx];
            u32 s2 = sol2[idx];
            uint8_t lp = cfg.log_primes[idx];
            if (s1 < B_size)
            {
                s_ptr[s1] += lp;
                in_block_large1.push_back(idx);
            }
            if (s2 < B_size)
            {
                s_ptr[s2] += lp;
                in_block_large2.push_back(idx);
            }
        }

        // 3. Поиск выживших кандидатов
        std::unordered_map<u32, std::vector<u32>> survivor_hits;
        for (u32 s_idx = 0; s_idx < B_size; ++s_idx)
        {
            if (s_ptr[s_idx] >= target_log)
            {
                survivor_hits[s_idx] = {};
            }
        }
        if (survivor_hits.empty())
            continue;

        for (u32 idx : in_block_large1)
        {
            auto it = survivor_hits.find(sol1[idx]);
            if (it != survivor_hits.end())
                it->second.push_back(idx);
        }
        for (u32 idx : in_block_large2)
        {
            auto it = survivor_hits.find(sol2[idx]);
            if (it != survivor_hits.end())
                it->second.push_back(idx);
        }

        // 4. Точное разложение выживших кандидатов
        for (const auto &kv : survivor_hits)
        {
            u32 s_idx = kv.first;
            const auto &large_divs = kv.second;

            // Применяем XOR-маску вокруг центра сита half_B
            u32 mask = half_B ^ s_idx;
            i64 y = (i64)(half_B ^ mask) - (i64)half_B;

            // X = A*y + B: считаем |X| и знак X с учётом знака B.
            bool y_neg = (y < 0);
            U512 ay = mul_u64(A, (u64)(y_neg ? -y : y));
            U512 X_raw;
            bool x_neg;
            if (y_neg == B_neg)
            {
                X_raw = add512(ay, B_mag);
                x_neg = y_neg;
            }
            else if (ay >= B_mag)
            {
                X_raw = sub512(ay, B_mag);
                x_neg = y_neg;
            }
            else
            {
                X_raw = sub512(B_mag, ay);
                x_neg = B_neg;
            }

            U512 X2 = mul512(X_raw, X_raw);
            if (X2 == cfg.kM)
                continue;

            int sign = 1;
            U512 diff;
            if (X2 >= cfg.kM)
            {
                diff = sub512(X2, cfg.kM);
                sign = 1;
            }
            else
            {
                diff = sub512(cfg.kM, X2);
                sign = -1;
            }

            U512 temp = div512(diff, A);
            std::vector<FactorExp> exps;
            exps.reserve(24);
            for (u32 a_idx : chosen_indices)
                exps.push_back({a_idx, 1});

            int tz = (temp.w[Z0] != 0) ? __builtin_ctzll(temp.w[Z0]) : 0;
            if ((temp.w[Z0] & 1) == 0)
            {
                u32 cnt2 = 0;
                while (!temp.is_zero() && (temp.w[Z0] & 1) == 0)
                {
                    int step = (temp.w[Z0] != 0) ? __builtin_ctzll(temp.w[Z0]) : 63;
                    temp.shr_k(step);
                    cnt2 += step;
                }
                exps.push_back({Z0, cnt2});
            }

            for (u32 idx : large_divs)
            {
                u32 p = cfg.prime_base[idx];
                if (temp.mod_u64(p) == 0)
                {
                    u32 cnt = 0;
                    while (temp.mod_u64(p) == 0)
                    {
                        temp.divrem_u64(p);
                        ++cnt;
                    }
                    merge_exps(exps, {{idx, cnt}});
                }
            }

            for (u32 idx = 1; idx < cfg.split_idx; ++idx)
            {
                u32 p = cfg.prime_base[idx];
                u32 rem_s = s_idx % p;
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
                        merge_exps(exps, {{idx, cnt}});
                        if (temp.bit_length() <= 1)
                            break;
                    }
                }
            }

            if (temp.bit_length() > 63)
                continue;
            u64 cof = temp.w[Z0];
            if (cof > cfg.lp2_bound)
                continue;

            if (cof > 1 && cof <= cfg.last_prime)
            {
                auto pit = cfg.prime_to_idx.find(cof);
                if (pit != cfg.prime_to_idx.end())
                {
                    merge_exps(exps, {{pit->second, 1}});
                    cof = 1;
                }
                else
                    continue;
            }

            U512 X_val = mod512(X_raw, cfg.M);
            if (x_neg && !X_val.is_zero()) {
                X_val = sub512(cfg.M, X_val);
            }

            if (cof == 1)
            {
                found_relations.push_back({1, 1, X_val, sign, std::move(exps)});
            }
            else if (cof <= cfg.lp_bound)
            {
                found_relations.push_back({1, cof, X_val, sign, std::move(exps)});
            }
            else if (cof >= cfg.min_2lp_comp)
            {
                if (!is_prime64(cof))
                {
                    u64 f1 = 0, f2 = 0;
                    if (brent_factor_two_primes(cof, f1, f2))
                    {
                        if (f1 > cfg.last_prime && f1 <= cfg.lp_bound &&
                            f2 > cfg.last_prime && f2 <= cfg.lp_bound)
                        {
                            if (is_prime64(f1) && is_prime64(f2))
                            {
                                found_relations.push_back({f1, f2, X_val, sign, std::move(exps)});
                            }
                        }
                    }
                }
            }
        }
    }
    return found_relations;
}

// Блок 6 из 6: Координатор `OpenMP` и функция `main` (205, 232 и 245 бит!)
// Здесь все 20 потоков вашего процессора параллельно запускают `process_one_A_batch_cpp` в цикле `#pragma omp parallel for schedule(dynamic, 1)`, без сериализации `pickle` и межпроцессных труб!

std::vector<u32> sieve_primes_cpp(u32 limit)
{
    std::vector<uint8_t> is_p(limit + 1, 1);
    is_p[Z0] = 0;
    if (limit >= 1)
        is_p[1] = 0;
    for (u32 i = 2; i * i <= limit; ++i)
    {
        if (is_p[i])
        {
            for (u32 j = i * i; j <= limit; j += i)
                is_p[j] = 0;
        }
    }
    std::vector<u32> primes;
    for (u32 i = 2; i <= limit; ++i)
        if (is_p[i])
            primes.push_back(i);
    return primes;
}

u32 find_best_multiplier_cpp(const U512 &M, const std::vector<u32> &small_primes)
{
    u32 best_k = 1;
    double best_score = -1e9;
    static const u32 cands[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 41, 42, 43, 46, 47};
    for (u32 k : cands)
    {
        U512 kM = mul_u64(M, k);
        double score = -0.5 * std::log((double)k);
        u64 mod8 = kM.w[Z0] & 7ULL;
        if (mod8 == 1)
            score += 2.0 * std::log(2.0);
        else if (mod8 == 5)
            score += 1.0 * std::log(2.0);
        for (size_t i = 1; i < std::min((size_t)40, small_primes.size()); ++i)
        {
            u32 p = small_primes[i];
            u64 rem = kM.mod_u64(p);
            if (rem == 0)
                score += std::log((double)p) / p;
            else if (pow_mod64(rem, (p - 1) / 2, p) == 1)
                score += 2.0 * std::log((double)p) / (p - 1);
        }
        if (score > best_score)
        {
            best_score = score;
            best_k = k;
        }
    }
    return best_k;
}

struct AutoParams {
    u32 prime_bound;
    int block_bits;
    int num_a_factors;
};

AutoParams auto_tune_siqs(const U512& M) {
    int bits = M.bit_length();
    double ln_M = bits * 0.6931471805599453; // bits * ln(2)
    double ln_ln_M = std::log(ln_M);
    
    // 1. Аппроксимация L_M[1/2, c] со степенной поправкой под 2LP SIQS.
    // Коэффициенты откалиброваны по эталонным прогонам:
    // 205 бит -> ~380k, 232 бита -> ~600k, 245 бит -> ~1.1M, 262 бита -> ~1.6M.
    double exponent = 0.357 * std::sqrt(ln_M * ln_ln_M) + 3.3;
    double raw_bound = std::exp(exponent);
    
    // Округляем до удобных десятков тысяч
    u32 prime_bound = (u32)(std::round(raw_bound / 10000.0) * 10000.0);
    prime_bound = std::max(50000U, std::min(prime_bound, 8000000U));

    // 2. Размер блока: держим в пределах L1/L2 кэша процессора (64-128 КБ)
    int block_bits = 16;
    if (bits >= 200) block_bits = 17; // 128 КБ массив (отлично помещается в L2)
    if (bits >= 270) block_bits = 18; // 256 КБ массив

    // 3. Расчет количества простых множителей s для коэффициента A
    // target_A ≈ sqrt(2 * M) / (B_size / 2)
    double half_B = (double)(1U << (block_bits - 1));
    double log2_target_A = 0.5 * (bits + 1.0) - std::log2(half_B);
    
    // Среднее простое из пула A выбирается около p ≈ prime_bound^0.65
    double log2_avg_q = std::log2((double)prime_bound) * 0.62;
    int s = (int)std::round(log2_target_A / log2_avg_q);
    s = std::max(5, std::min(s, 10));

    return {prime_bound, block_bits, s};
}


void factor_siqs_cpp(const std::string &dec_M, u32 prime_bound, int block_bits, int num_a_factors)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    U512 M = from_dec_string(dec_M);
    auto raw_primes = sieve_primes_cpp(prime_bound);
    u32 k = find_best_multiplier_cpp(M, raw_primes);
    U512 kM = mul_u64(M, k);

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
            u32 r = (u32)(kM.w[Z0] & 1ULL);
            cfg.root1.push_back(r);
            cfg.root2.push_back(r);
            cfg.log_primes.push_back(1);
        }
        else
        {
            u64 rem = kM.mod_u64(p);
            if (rem == 0 || pow_mod64(rem, (p - 1) / 2, p) == 1)
            {
                u32 r = tonelli_shanks(rem, p);
                cfg.prime_base.push_back(p);
                cfg.root1.push_back(r);
                cfg.root2.push_back((r != p - r) ? (p - r) : r);
                cfg.log_primes.push_back((uint8_t)std::llround(std::log2((double)p)));
            }
        }
    }

    cfg.num_primes = (u32)cfg.prime_base.size();
    for (u32 i = 0; i < cfg.num_primes; ++i)
        cfg.prime_to_idx[cfg.prime_base[i]] = i;
    cfg.last_prime = cfg.prime_base.back();
    cfg.min_2lp_comp = cfg.last_prime * cfg.last_prime;
    cfg.lp_bound = (u64)prime_bound * 30ULL;
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

    U512 two_kM = kM;
    two_kM.shl1();
    U512 target_A = div512(isqrt512(two_kM), U512(cfg.half_B));
    double log2_A = (double)target_A.bit_length();
    u32 ideal_q = (u32)std::max(50.0, std::pow(2.0, log2_A / cfg.s));

    u32 center_idx = cfg.num_primes / 2;
    for (u32 i = 0; i < cfg.num_primes; ++i)
    {
        if (cfg.prime_base[i] >= ideal_q)
        {
            center_idx = i;
            break;
        }
    }
    u32 pool_start = std::max(15U, std::min(center_idx > 60 ? center_idx - 60 : 15U, cfg.num_primes - 130));
    for (u32 i = pool_start; i < std::min(cfg.num_primes - 1, pool_start + 120); ++i)
    {
        cfg.pool_indices.push_back(i);
    }

    int threads = omp_get_max_threads();
    std::cout << "\n=== C++ NATIVE 2LP SIQS (" << M.bit_length() << " бит, потоков OpenMP: " << threads << ") ===\n";
    std::cout << "M = " << dec_M << "\n";
    std::cout << "Множитель k = " << k << ", База: " << cfg.num_primes << " простых (до " << prime_bound << ")\n";

    MarkowitzGaussianGF2 gauss;
    LargePrimeGraph lp_graph(M);
    u64 total_polys = 0;
    u64 batch_seed = 1000;
    u64 polys_per_A = 1ULL << (cfg.s - 1);

    while (true)
    {
        std::vector<std::vector<RawRelation>> batch_results(threads);

#pragma omp parallel for schedule(dynamic, 1)
        for (int t = 0; t < threads; ++t)
        {
            batch_results[t] = process_one_A_batch_cpp(cfg, batch_seed + t + 1);
        }
        batch_seed += threads;
        total_polys += threads * polys_per_A;

        for (int t = 0; t < threads; ++t)
        {
            for (const auto &r : batch_results[t])
            {
                FullRelation full_rel;
                if (lp_graph.add_edge_and_extract_cycle(r, full_rel))
                {
                    gauss.add_relation(std::move(full_rel));
                }
            }
        }

        auto now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(now - t0).count();
        std::cout << "[" << elapsed << " с | Полином " << total_polys << "] Соотношений: " << gauss.relations.size()
                  << " | Простых (>=2 вх.): " << gauss.multi_hit_primes
                  << " (синглетонов: " << gauss.singleton_primes << ", циклов 2LP: " << lp_graph.cycles_found << ")\n";

        if (gauss.ready_to_check_rank())
        {
            auto tg0 = std::chrono::high_resolution_clock::now();
            auto deps = gauss.solve_pruned_matrix();
            for (const auto &dep : deps)
            {
                U512 LHS(1);
                std::unordered_map<u32, u32> total_exps;
                for (u32 idx_r : dep)
                {
                    LHS = mul_mod512(LHS, gauss.relations[idx_r].X_val, M);
                    for (const auto &fe : gauss.relations[idx_r].exps)
                    {
                        total_exps[fe.idx] += fe.cnt;
                    }
                }
                U512 RHS(1);
                for (const auto &kv : total_exps)
                {
                    u32 p_val = cfg.prime_base[kv.first];
                    U512 p_pow = pow_mod512(U512(p_val), kv.second / 2, M);
                    RHS = mul_mod512(RHS, p_pow, M);
                }
                U512 diff = (LHS >= RHS) ? sub512(LHS, RHS) : sub512(RHS, LHS);
                U512 sum = add512(LHS, RHS);
                if (sum >= M)
                    sum.sub(M);

                for (const U512 &cand : {diff, sum})
                {
                    U512 g = gcd512(cand, M);
                    if (g > U512(1) && g < M)
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        double tot_s = std::chrono::duration<double>(t1 - t0).count();
                        double g_s = std::chrono::duration<double>(t1 - tg0).count();
                        U512 q_fac = div512(M, g);
                        std::cout << "[УСПЕХ C++ за " << tot_s << " сек (Гаусс: " << g_s << " сек)!]\n";
                        std::cout << "Полиномов: " << total_polys << ", Проверено XOR-масок: " << (total_polys * cfg.B_size) << "\n";
                        std::cout << "Найденные делители:\n  p = " << to_dec_string(g) << "\n  q = " << to_dec_string(q_fac) << "\n";
                        return;
                    }
                }
            }
        }
    }
}

void factor_siqs_auto(const std::string& dec_M) {
    U512 M = from_dec_string(dec_M);
    AutoParams p = auto_tune_siqs(M);
    std::cout << "[Автоподбор] Бит: " << M.bit_length()
              << " -> prime_bound = " << p.prime_bound
              << ", block_bits = " << p.block_bits
              << ", num_a_factors = " << p.num_a_factors << "\n";
    
    factor_siqs_cpp(dec_M, p.prime_bound, p.block_bits, p.num_a_factors);
}



int main()
{
    // 205-битное число (в Python было 22.15 сек)
    // factor_siqs_cpp("33801295533024112192430315909773114733482083162238224697211891", 380000, 17, 7);

    // 232-битное число (в Python было 98.6 сек)
    // factor_siqs_cpp("5911987939788655955220092655904964287465054776969341694771422110974013", 600000, 17, 8);

    // 245-битное число (74 знака. В Python было 272.2 сек):
    // 5183690555951629707790497436756762363 * 5864810964374268616720755200996315407
    // factor_siqs_cpp("30401364784979048844780572496488416288730829724034373452010539392784830241", 1100000, 17, 8);

    std::string M_262 = "5995806622530429600530984275000092882631332762476331275490451238005952121568939";    
    factor_siqs_auto(M_262);

    return 0;
}
