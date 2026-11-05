#pragma once

// ---------------------------------------------------------------------------
// siqs/gauss.hpp — разреженный гаусс над GF(2).
//
// Исключение с XOR-историей строк: на выходе не вектор решения, а набор
// строк-зависимостей, по которым потом считается кандидат в делители.
// ---------------------------------------------------------------------------

#include "common.hpp"
#include "relation.hpp"

namespace bignum::siqs_xor {

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
        // Нужен положительный избыток строк над активными столбцами. Прежнее
        // условие (n_rels + 120 >= total_active) разрешало запуск при дефиците
        // до 120 столбцов: исключение заведомо не давало зависимости, а
        // last_check_count откладывал следующую попытку на 500 соотношений —
        // на 160 битах это ~0.25 с впустую.
        const size_t margin = std::max<size_t>(40, total_active / 50);
        return n_rels >= total_active + margin;
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

} // namespace bignum::siqs_xor
