import random
import math
import time
import multiprocessing as mp

# --- 1. Решето Эратосфена и генерация простых ---
def sieve_primes(limit):
    ZERO, ONE = int(0), int(1)
    is_prime = list(range(limit + 1))
    is_prime[ZERO] = 0
    if limit >= 1:
        is_prime[ONE] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i] != 0:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = 0
    return [p for p in is_prime if p != 0]

def is_probable_prime(n, k=4):
    if n < 2: return False
    if n in (2, 3, 5, 7): return True
    if (n & 1) == 0 or n % 3 == 0: return False
    r, d = 0, n - 1
    while (d & 1) == 0:
        r += 1
        d >>= 1
    for a in (2, 3, 5, 7, 11, 13)[:k]:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True

def generate_prime(bits):
    while True:
        p = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_probable_prime(p, k=6):
            return p

# --- 2. Быстрый метод Брента для расщепления 2LP ---
def brent_factor_two_primes(n):
    if (n & 1) == 0:
        return 2, n >> 1
    for c in (1, 3, 5):
        y = 2
        m = 32
        g = 1
        r = 1
        q = 1
        ys = 2
        x = 2
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                lim = min(m, r - k)
                for _ in range(lim):
                    y = (y * y + c) % n
                    q = (q * abs(x - y)) % n
                g = math.gcd(q, n)
                k += m
            r <<= 1
            if r > 512:
                break
        if g == n:
            while True:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
                if g > 1:
                    break
        if 1 < g < n:
            return g, n // g
    return None, None

# --- 3. Алгоритм Тонелли — Шэнкса и множитель Кнута — Шрёппеля ---
def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if p == 2: return n
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while (q & 1) == 0:
        q >>= 1
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, temp = 1, (t * t) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, (b * b) % p, (t * b * b) % p, (r * b) % p
    return r

def find_best_multiplier(M, small_primes):
    best_k = 1
    best_score = -1e9
    candidates = (1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 41, 42, 43, 46, 47)
    for k in candidates:
        kM = k * M
        score = -0.5 * math.log(k)
        if kM % 8 == 1:
            score += 2.0 * math.log(2)
        elif kM % 8 == 5:
            score += 1.0 * math.log(2)
        for p in small_primes[1:40]:
            if kM % p == 0:
                score += math.log(p) / p
            elif pow(kM % p, (p - 1) // 2, p) == 1:
                score += 2.0 * math.log(p) / (p - 1)
        if score > best_score:
            best_score = score
            best_k = k
    return best_k

# --- 4. Решатель Гаусса Марковица (Глубокое сжатие весов 2..5 + Быстрое ядро) ---
class MarkowitzGaussianGF2:
    def __init__(self, num_primes):
        self.num_primes = num_primes
        self.relations = []
        self.col_freq = {}
        self.multi_hit_primes = 0
        self.singleton_primes = 0
        self.last_check_count = 0

    def add_relation(self, X_val, sign, exps):
        self.relations.append((X_val, sign, exps))
        for p_idx, cnt in exps.items():
            if cnt & 1:
                prev = self.col_freq.get(p_idx, 0)
                if prev == 0:
                    self.singleton_primes += 1
                elif prev == 1:
                    self.singleton_primes -= 1
                    self.multi_hit_primes += 1
                self.col_freq[p_idx] = prev + 1

    def ready_to_check_rank(self):
        n_rels = len(self.relations)
        if n_rels - self.last_check_count < 500:
            return False
        total_active = self.multi_hit_primes + self.singleton_primes
        return n_rels >= total_active - 120

    def solve_pruned_matrix(self):
        self.last_check_count = len(self.relations)
        ZERO, ONE = int(0), int(1)

        active_rows = {}
        col_to_rows = {}
        for r_idx, (_, sign, exps) in enumerate(self.relations):
            odd_cols = {p_idx for p_idx, cnt in exps.items() if (cnt & 1)}
            s_bit = 1 if sign == -1 else 0
            active_rows[r_idx] = [odd_cols, s_bit, 1 << r_idx]
            for c in odd_cols:
                if c not in col_to_rows:
                    col_to_rows[c] = {r_idx}
                else:
                    col_to_rows[c].add(r_idx)

        def remove_row(r_del):
            cols_del = active_rows[r_del][ZERO]
            del active_rows[r_del]
            for c in cols_del:
                s_set = col_to_rows.get(c)
                if s_set is not None:
                    s_set.discard(r_del)
                    if not s_set:
                        del col_to_rows[c]

        def xor_row_into(r_src, r_dst):
            src_cols, src_s, src_h = active_rows[r_src]
            dst_entry = active_rows[r_dst]
            dst_cols = dst_entry[ZERO]
            for c in src_cols:
                if c in dst_cols:
                    dst_cols.remove(c)
                    col_to_rows[c].discard(r_dst)
                    if not col_to_rows[c]:
                        del col_to_rows[c]
                else:
                    dst_cols.add(c)
                    if c not in col_to_rows:
                        col_to_rows[c] = {r_dst}
                    else:
                        col_to_rows[c].add(r_dst)
            dst_entry[ONE] ^= src_s
            dst_entry[int(2)] ^= src_h

        def prune_singletons():
            singleton_queue = [c for c, r_set in col_to_rows.items() if len(r_set) == 1]
            while singleton_queue:
                c = singleton_queue.pop()
                r_set = col_to_rows.get(c)
                if r_set and len(r_set) == 1:
                    r_del = next(iter(r_set))
                    cols_affected = list(active_rows[r_del][ZERO])
                    remove_row(r_del)
                    for c_aff in cols_affected:
                        s_aff = col_to_rows.get(c_aff)
                        if s_aff and len(s_aff) == 1:
                            singleton_queue.append(c_aff)

        # 1. Начальная очистка синглетонов
        prune_singletons()
        if len(active_rows) <= len(col_to_rows) + 5:
            return []

        orig_core_rows = len(active_rows)

        # 2. Многопроходное каскадное схлопывание столбцов веса 2, 3, 4 и 5 (Structured Gaussian Elimination)
        for max_weight in (2, 2, 3, 3, 4, 5):
            candidates = [c for c, r_set in col_to_rows.items() if 2 <= len(r_set) <= max_weight]
            # Сначала схлопываем столбцы с самыми короткими строками, чтобы не раздувать вес!
            for c in candidates:
                r_set = col_to_rows.get(c)
                if not r_set or len(r_set) > max_weight or len(r_set) < 2:
                    continue
                if len(active_rows) <= len(col_to_rows) + 24:
                    break
                r_list = sorted(r_set, key=lambda r: len(active_rows[r][ZERO]))
                pivot_r = r_list[ZERO]
                # Ограничиваем вес опорной строки, чтобы избежать взрыва заполнения (fill-in)
                if len(active_rows[pivot_r][ZERO]) > 35:
                    continue
                for target_r in r_list[ONE:]:
                    xor_row_into(pivot_r, target_r)
                remove_row(pivot_r)
            prune_singletons()

        print(f"    [Гаусс-триггер] Ядро после синглетонов: {orig_core_rows} строк -> "
              f"сжато SGE (веса 2..5) до {len(active_rows)} строк и {len(col_to_rows)} столбцов!")

        # 3. Быстрый битовый Гаусс Марковица на компактном ядре
        sorted_cols = sorted(col_to_rows.keys(), key=lambda c: len(col_to_rows[c]), reverse=True)
        col_map = {c: i for i, c in enumerate(sorted_cols)}
        sign_bit = len(col_map)

        sorted_row_indices = sorted(active_rows.keys(), key=lambda r: len(active_rows[r][ZERO]))
        pivot_map = {}
        dependencies = []

        for r_idx in sorted_row_indices:
            cols, s_bit, history = active_rows[r_idx]
            vec = (1 << sign_bit) if s_bit else 0
            for c in cols:
                if c in col_map:
                    vec |= (1 << col_map[c])

            while vec > 0:
                msb = vec.bit_length() - 1
                if msb in pivot_map:
                    p_vec, p_hist = pivot_map[msb]
                    vec ^= p_vec
                    history ^= p_hist
                else:
                    pivot_map[msb] = (vec, history)
                    break

            if vec == 0 and history > 0:
                dep = []
                temp_h = history
                while temp_h:
                    lowbit = temp_h & -temp_h
                    dep.append(lowbit.bit_length() - 1)
                    temp_h ^= lowbit
                dependencies.append(dep)
                if len(dependencies) >= 16:
                    break

        return dependencies

# --- 5. Сверхбыстрый граф циклов 2LP на ориентированном остовном лесе ---
class LargePrimeGraph:
    def __init__(self, M):
        self.M = M
        self.uf_parent = {}
        self.tree_parent = {}
        self.num_edges = 0
        self.cycles_found = 0

    def find_root(self, v):
        curr = v
        while curr in self.uf_parent and self.uf_parent[curr] != curr:
            curr = self.uf_parent[curr]
        root = curr
        curr = v
        while curr in self.uf_parent and self.uf_parent[curr] != root:
            nxt = self.uf_parent[curr]
            self.uf_parent[curr] = root
            curr = nxt
        if root not in self.uf_parent:
            self.uf_parent[root] = root
        return root

    def _reroot_tree(self, node):
        path = []
        curr = node
        while curr in self.tree_parent:
            path.append((curr, self.tree_parent[curr]))
            curr = self.tree_parent[curr][int(0)]
        for child, (par, e_X, e_sign, e_exps, e_primes) in reversed(path):
            self.tree_parent[par] = (child, e_X, e_sign, e_exps, e_primes)
        if node in self.tree_parent:
            del self.tree_parent[node]

    def add_edge_and_extract_cycle(self, u, v, X_val, sign, exps):
        if u == v:
            inv_u = pow(u, -1, self.M)
            return (X_val * inv_u) % self.M, sign, exps

        ru = self.find_root(u)
        rv = self.find_root(v)
        edge_primes = [p for p in (u, v) if p > 1]

        if ru != rv:
            self._reroot_tree(u)
            self.tree_parent[u] = (v, X_val, sign, exps, edge_primes)
            self.uf_parent[ru] = rv
            self.num_edges += 1
            return None

        anc_u = {}
        curr = u
        while True:
            anc_u[curr] = True
            if curr not in self.tree_parent:
                break
            curr = self.tree_parent[curr][int(0)]

        curr = v
        while curr not in anc_u:
            curr = self.tree_parent[curr][int(0)]
        lca = curr

        comb_X = X_val
        comb_sign = sign
        comb_exps = exps.copy()
        lp_counts = {}
        for p in edge_primes:
            lp_counts[p] = lp_counts.get(p, 0) + 1

        for start_node in (u, v):
            curr = start_node
            while curr != lca:
                par, e_X, e_sign, e_exps, e_primes = self.tree_parent[curr]
                comb_X = (comb_X * e_X) % self.M
                comb_sign *= e_sign
                for k_p, v_p in e_exps.items():
                    comb_exps[k_p] = comb_exps.get(k_p, 0) + v_p
                for p in e_primes:
                    lp_counts[p] = lp_counts.get(p, 0) + 1
                curr = par

        denom = 1
        for lp_val, cnt in lp_counts.items():
            denom = (denom * pow(lp_val, cnt // 2, self.M)) % self.M
        inv_denom = pow(denom, -1, self.M)
        comb_X = (comb_X * inv_denom) % self.M

        self.cycles_found += 1
        return comb_X, comb_sign, comb_exps

# --- 6. Воркер SIQS с Active-Block Bucket (в 5 раз меньше проверок крупных простых со 100% точностью!) ---
GLOBAL_SIQS_DATA = {}

def init_siqs_worker(data_dict):
    global GLOBAL_SIQS_DATA
    GLOBAL_SIQS_DATA = data_dict

def process_one_A_batch_2lp(seed_offset):
    ZERO, ONE = int(0), int(1)
    d = GLOBAL_SIQS_DATA
    M = d['M']
    kM = d['kM']
    prime_base = d['prime_base']
    prime_to_idx = d['prime_to_idx']
    roots = d['roots']
    log_primes = d['log_primes']
    num_primes = d['num_primes']
    last_prime = d['last_prime']
    min_2lp_comp = last_prime * last_prime
    lp_bound = d['lp_bound']
    lp2_bound = d['lp2_bound']
    B_size = d['B_size']
    half_B = d['half_B']
    s = d['s']
    pool_indices = d['pool_indices']
    small_prime_skip = d['small_prime_skip']
    split_idx = d['split_idx']
    trans_tables = d['trans_tables']
    zero_chunk = b'\x00' * 512

    rng = random.Random(seed_offset ^ int(time.time() * 1000000))
    chosen_indices = sorted(rng.sample(pool_indices, s))
    chosen_set = set(chosen_indices)
    A = 1
    for idx in chosen_indices:
        A *= prime_base[idx]

    B_terms = []
    for idx in chosen_indices:
        q_j = prime_base[idx]
        r_j = roots[idx][ZERO]
        A_div_q = A // q_j
        inv = pow(A_div_q % q_j, -1, q_j)
        gamma = (r_j * inv) % q_j
        if gamma > q_j // 2:
            gamma = q_j - gamma
        B_terms.append(A_div_q * gamma)

    B = sum(B_terms)
    inv_A_mod_p = [ZERO] * num_primes
    B_terms_mod_p = [[ZERO] * num_primes for _ in range(s)]
    sol1 = [ZERO] * num_primes
    sol2 = [ZERO] * num_primes

    for idx in range(1, num_primes):
        if idx in chosen_set:
            continue
        p = prime_base[idx]
        inv_A = pow(A % p, -1, p)
        inv_A_mod_p[idx] = inv_A
        for j in range(s):
            B_terms_mod_p[j][idx] = (2 * B_terms[j] * inv_A) % p

        r_list = roots[idx]
        r1 = r_list[ZERO]
        r2 = r_list[ONE] if len(r_list) > 1 else r1
        B_mod = B % p
        sol1[idx] = ((r1 - B_mod) * inv_A + half_B) % p
        sol2[idx] = ((r2 - B_mod) * inv_A + half_B) % p

    signs = [ONE] * s
    polys_per_A = 1 << (s - 1)
    found_relations = []

    for poly_i in range(polys_per_A):
        if poly_i > 0:
            gray_curr = poly_i ^ (poly_i >> 1)
            gray_prev = (poly_i - 1) ^ ((poly_i - 1) >> 1)
            diff = gray_curr ^ gray_prev
            v = diff.bit_length()

            new_sign = -signs[v]
            signs[v] = new_sign
            row_v = B_terms_mod_p[v]

            if new_sign == -1:
                B -= 2 * B_terms[v]
                for idx in range(1, num_primes):
                    if inv_A_mod_p[idx]:
                        p = prime_base[idx]
                        delta = row_v[idx]
                        s1_new = sol1[idx] + delta
                        sol1[idx] = s1_new - p if s1_new >= p else s1_new
                        s2_new = sol2[idx] + delta
                        sol2[idx] = s2_new - p if s2_new >= p else s2_new
            else:
                B += 2 * B_terms[v]
                for idx in range(1, num_primes):
                    if inv_A_mod_p[idx]:
                        p = prime_base[idx]
                        delta = row_v[idx]
                        s1_new = sol1[idx] - delta
                        sol1[idx] = s1_new + p if s1_new < 0 else s1_new
                        s2_new = sol2[idx] - delta
                        sol2[idx] = s2_new + p if s2_new < 0 else s2_new

        C = (B * B - kM) // A
        sieve_logs = bytearray(B_size)

        # 1. Просеиваем малые и средние простые (p < B_size)
        for idx in range(small_prime_skip, split_idx):
            if inv_A_mod_p[idx] == 0:
                continue
            p = prime_base[idx]
            lp = log_primes[idx]
            s1 = sol1[idx]
            s2 = sol2[idx]
            for s_idx in range(s1, B_size, p):
                sieve_logs[s_idx] += lp
            if s2 != s1:
                for s_idx in range(s2, B_size, p):
                    sieve_logs[s_idx] += lp

        # 2. Просеиваем крупные простые (p >= B_size) и запоминаем ТОЛЬКО те, что реально попали в блок < B_size!
        in_block_large1 = []
        in_block_large2 = []
        for idx in range(split_idx, num_primes):
            s1 = sol1[idx]
            s2 = sol2[idx]
            lp = log_primes[idx]
            if s1 < B_size:
                sieve_logs[s1] += lp
                in_block_large1.append(idx)
            if s2 < B_size:
                sieve_logs[s2] += lp
                in_block_large2.append(idx)

        sample_Fy = abs(C) if C != 0 else (half_B * math.isqrt(2 * kM))
        target_log = min(250, max(15, sample_Fy.bit_length() - int(math.log2(lp2_bound)) + 3))

        filtered = sieve_logs.translate(trans_tables[target_log])

        # 3. Собираем выжившие кандидаты блоками по 512 байт
        survivor_hits = {}
        for chunk_start in range(0, B_size, 512):
            chunk = filtered[chunk_start : chunk_start + 512]
            if chunk == zero_chunk:
                continue
            for offset, s_log_val in enumerate(chunk):
                if s_log_val != 0:
                    survivor_hits[chunk_start + offset] = []

        if not survivor_hits:
            continue

        # 4. Проверяем ТОЛЬКО активные крупные простые из in_block_large (в 5-6 раз меньше итераций!)
        for idx in in_block_large1:
            s1 = sol1[idx]
            if s1 in survivor_hits:
                survivor_hits[s1].append(idx)
        for idx in in_block_large2:
            s2 = sol2[idx]
            if s2 in survivor_hits:
                survivor_hits[s2].append(idx)

        # 5. Разложение выживших кандидатов
        for s_idx, large_div_indices in survivor_hits.items():
            mask = half_B ^ s_idx
            y = (half_B ^ mask) - half_B
            X_raw = A * y + B
            Fy = (X_raw * X_raw - kM) // A
            if Fy == 0:
                continue

            sign = -1 if Fy < 0 else 1
            temp = abs(Fy)
            exps = {}
            for a_idx in chosen_indices:
                exps[a_idx] = 1

            if (temp & 1) == 0:
                cnt = 0
                while (temp & 1) == 0:
                    temp >>= 1
                    cnt += 1
                exps[ZERO] = cnt

            for idx in large_div_indices:
                p = prime_base[idx]
                if temp % p == 0:
                    cnt = 0
                    while temp % p == 0:
                        temp //= p
                        cnt += 1
                    exps[idx] = exps.get(idx, 0) + cnt

            for idx in range(1, split_idx):
                p = prime_base[idx]
                rem_s = s_idx % p
                if rem_s == sol1[idx] or rem_s == sol2[idx]:
                    if temp % p == 0:
                        cnt = 0
                        while temp % p == 0:
                            temp //= p
                            cnt += 1
                        exps[idx] = exps.get(idx, 0) + cnt
                        if temp == 1:
                            break
                elif idx in chosen_set:
                    if temp % p == 0:
                        cnt = 0
                        while temp % p == 0:
                            temp //= p
                            cnt += 1
                        exps[idx] = exps.get(idx, 0) + cnt

            if temp > lp2_bound:
                continue

            if 1 < temp <= last_prime:
                found_idx = prime_to_idx.get(temp, -1)
                if found_idx != -1:
                    exps[found_idx] = exps.get(found_idx, 0) + 1
                    temp = 1
                else:
                    continue

            X_val = X_raw % M

            if temp == 1:
                found_relations.append((1, 1, X_val, sign, exps))
            elif temp <= lp_bound:
                found_relations.append((1, temp, X_val, sign, exps))
            elif temp >= min_2lp_comp:
                if not is_probable_prime(temp, k=2):
                    f1, f2 = brent_factor_two_primes(temp)
                    if f1 is not None and last_prime < f1 <= lp_bound and last_prime < f2 <= lp_bound:
                        if is_probable_prime(f1, k=2) and is_probable_prime(f2, k=2):
                            found_relations.append((f1, f2, X_val, sign, exps))

    return polys_per_A, found_relations

# --- 7. Главный координатор: Hybrid Hit-Map + Early-Rank Markowitz Gauss ---
def factor_siqs_245bit_edition(M, prime_bound=850000, block_bits=17, num_a_factors=8, workers=None):
    t0 = time.time()
    ZERO, ONE = int(0), int(1)
    raw_primes = sieve_primes(prime_bound)
    k = find_best_multiplier(M, raw_primes)
    kM = k * M

    prime_base = []
    roots = []
    log_primes = []
    for p in raw_primes:
        if p == 2:
            prime_base.append(2)
            roots.append((kM % 2,))
            log_primes.append(1)
        else:
            rem = kM % p
            if rem == 0 or pow(rem, (p - 1) // 2, p) == 1:
                r = tonelli_shanks(rem, p)
                prime_base.append(p)
                roots.append((r, p - r) if r != p - r else (r,))
                log_primes.append(int(round(math.log2(p))))

    num_primes = len(prime_base)
    prime_to_idx = {p: i for i, p in enumerate(prime_base)}
    last_prime = prime_base[num_primes - 1]
    lp_bound = prime_bound * 30
    lp2_bound = lp_bound * lp_bound

    B_size = 1 << block_bits
    half_B = B_size // 2

    split_idx = num_primes
    for idx, p in enumerate(prime_base):
        if p >= B_size:
            split_idx = idx
            break

    target_A = math.isqrt(2 * kM) // half_B
    s = num_a_factors
    ideal_q = max(50, int(target_A ** (1.0 / s)))

    center_idx = num_primes // 2
    for idx, p in enumerate(prime_base):
        if p >= ideal_q:
            center_idx = idx
            break
    pool_start = max(15, min(center_idx - 60, num_primes - 130))
    pool_indices = list(range(pool_start, min(num_primes - 1, pool_start + 120)))

    trans_tables = {}
    for th in range(1, 255):
        tbl = bytearray(256)
        for b_val in range(th, 256):
            tbl[b_val] = b_val
        trans_tables[th] = bytes(tbl)

    if workers is None:
        workers = max(1, mp.cpu_count())

    shared_data = {
        'M': M, 'kM': kM,
        'prime_base': prime_base, 'prime_to_idx': prime_to_idx,
        'roots': roots, 'log_primes': log_primes,
        'num_primes': num_primes, 'last_prime': last_prime,
        'lp_bound': lp_bound, 'lp2_bound': lp2_bound,
        'B_size': B_size, 'half_B': half_B, 's': s,
        'pool_indices': pool_indices,
        'small_prime_skip': 17,
        'split_idx': split_idx,
        'trans_tables': trans_tables
    }

    print(f"=== HYBRID HIT-MAP + MARKOWITZ-GAUSS SIQS ({M.bit_length()} бит, ядер: {workers}) ===")
    print(f"M = {M}")
    print(f"Множитель k = {k}, База: {num_primes} простых (до {prime_bound}, малых до B_size: {split_idx})")
    print(f"Граница 1LP: {lp_bound:.2e}, Граница 2LP (L1*L2): {lp2_bound:.2e}")

    gauss = MarkowitzGaussianGF2(num_primes)
    lp_graph = LargePrimeGraph(M)
    total_polys = 0
    last_print_polys = 0

    pool = mp.Pool(processes=workers, initializer=init_siqs_worker, initargs=(shared_data,))

    try:
        seed_stream = range(1, 5000000)
        for polys_done, rel_list in pool.imap_unordered(process_one_A_batch_2lp, seed_stream, chunksize=1):
            total_polys += polys_done
            for u_lp, v_lp, X_val, sign, exps in rel_list:
                if u_lp == 1 and v_lp == 1:
                    full_rel = (X_val, sign, exps)
                else:
                    for lp_cand in (u_lp, v_lp):
                        if lp_cand > 1:
                            g_lp = math.gcd(lp_cand, M)
                            if 1 < g_lp < M:
                                return g_lp, M // g_lp
                    full_rel = lp_graph.add_edge_and_extract_cycle(u_lp, v_lp, X_val, sign, exps)

                if full_rel is None:
                    continue

                r_X, r_sign, r_exps = full_rel
                gauss.add_relation(r_X, r_sign, r_exps)

            if total_polys - last_print_polys >= 1280:
                last_print_polys = total_polys
                elapsed = time.time() - t0
                print(f"[{elapsed:5.1f} с | Полином {total_polys:5d}] Соотношений: {len(gauss.relations):5d} | "
                      f"Простых (>=2 вх.): {gauss.multi_hit_primes:5d} (синглетонов: {gauss.singleton_primes}, циклов 2LP: {lp_graph.cycles_found})")

            if gauss.ready_to_check_rank():
                t_g = time.time()
                deps = gauss.solve_pruned_matrix()
                for dep_indices in deps:
                    LHS = 1
                    total_exps = {}
                    for idx_r in dep_indices:
                        x_i, _, exp_i = gauss.relations[idx_r]
                        LHS = (LHS * x_i) % M
                        for p_idx, cnt in exp_i.items():
                            total_exps[p_idx] = total_exps.get(p_idx, 0) + cnt
                    RHS = 1
                    for p_idx, cnt in total_exps.items():
                        p_val = prime_base[p_idx]
                        RHS = (RHS * pow(p_val, cnt // 2, M)) % M
                    for cand in (LHS - RHS, LHS + RHS):
                        g = math.gcd(cand, M)
                        if 1 < g < M:
                            elapsed = time.time() - t0
                            g_time = time.time() - t_g
                            total_masks = total_polys * B_size
                            print(f"[УСПЕХ за {elapsed:.3f} сек (Гаусс Марковица: {g_time:.3f} сек)!]")
                            print(f"Полиномов: {total_polys}, Проверено XOR-масок: {total_masks}")
                            print(f"Полных соотношений: {len(gauss.relations)} / Простых (>=2 вх.): {gauss.multi_hit_primes} "
                                  f"(синглетонов отброшено: {gauss.singleton_primes}, циклов 2LP: {lp_graph.cycles_found})")
                            print(f"Найденные делители:\n  p = {g}\n  q = {M // g}\n")
                            return g, M // g
    finally:
        pool.terminate()
        pool.join()

if __name__ == "__main__":
    mp.freeze_support()

    # 1. Ваше 232-битное число с оптимальной базой 600 000 (вместо 850 000) и без холостых вызовов Гаусса!
    M_232 = 5911987939788655955220092655904964287465054776969341694771422110974013
    factor_siqs_245bit_edition(M_232, prime_bound=600000, block_bits=17, num_a_factors=8)

    # 2. Штурм 245-битного полупростого числа (74 десятичных знака!)
    p_122 = generate_prime(122)
    q_123 = generate_prime(123)
    factor_siqs_245bit_edition(p_122 * q_123, prime_bound=1100000, block_bits=17, num_a_factors=8)
