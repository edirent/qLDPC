import numpy as np

def _gf_mult(a, b, modulus):
    res = 0
    while b > 0:
        if b & 1:
            res ^= a
        a <<= 1
        if a & modulus[1]:
            a ^= modulus[0]
        b >>= 1
    return res

def _gf_add(a, b):
    return a ^ b

def _gf_inverse(a, modulus, field_size):
    if a == 0:
        raise ValueError("Cannot compute inverse of 0")
    res = 1
    power = field_size - 2
    base = a
    while power > 0:
        if power & 1:
            res = _gf_mult(res, base, modulus)
        base = _gf_mult(base, base, modulus)
        power >>= 1
    return res

def _find_irreducible_polynomial(degree):
    irreducible_polys = {
        1: 0b11,
        2: 0b111,
        3: 0b1011,
        4: 0b10011,
        5: 0b100101,
        6: 0b1000011,
        7: 0b10010011,
        8: 0b110110001,
    }
    if degree not in irreducible_polys:
        if degree == 8:
            return (0b100011101, 1 << degree)
        raise ValueError(f"No irreducible polynomial for degree {degree}.")
    return (irreducible_polys[degree], 1 << degree)

def _get_ord_2(n):
    if n == 1: return 1
    k, res = 1, (2 % n)
    while res != 1:
        res = (res * 2) % n
        k += 1
    return k

def two_pow_factor(n):
    odd, twop = n, 1
    while odd % 2 == 0:
        odd //= 2
        twop *= 2
    return odd, twop

def min_ext_deg(l, m):
    l1, _ = two_pow_factor(l)
    m1, _ = two_pow_factor(m)
    if l1 == m1 == 1:
        return 1
    return int(np.lcm(_get_ord_2(l1), _get_ord_2(m1)))

def _power_in_gf(base, exp, modulus, field_size):
    if exp == 0:
        return 1
    res, cur = 1, base
    while exp > 0:
        if exp & 1:
            res = _gf_mult(res, cur, modulus)
        cur = _gf_mult(cur, cur, modulus)
        exp >>= 1
    return res

def k_bb_code_final(l, m, A_exp, B_exp):
    n = 2 * l * m
    l1, le = two_pow_factor(l)
    m1, me = two_pow_factor(m)
    r = min_ext_deg(l1, m1) or 1

    if r == 1:
        modulus = (0, 0)
        field_elems = [0, 1]
    else:
        modulus = _find_irreducible_polynomial(r)
        field_elems = list(range(1, 1 << r))

    alphas = [1] if l1 == 1 else [a for a in field_elems if _power_in_gf(a, l1, modulus, 1<<r) == 1]
    betas  = [1] if m1 == 1 else [b for b in field_elems if _power_in_gf(b, m1, modulus, 1<<r) == 1]

    total_rank = 0

    def rank_gf2r(matrix, modulus_poly, field_size):
        rows, cols = matrix.shape
        rank = 0
        M = np.copy(matrix)
        for j in range(cols):
            piv = rank
            while piv < rows and M[piv, j] == 0:
                piv += 1
            if piv == rows:
                continue
            if piv != rank:
                M[[rank, piv]] = M[[piv, rank]]
            pv = M[rank, j]
            inv = _gf_inverse(pv, modulus_poly, field_size)
            for c in range(j, cols):
                M[rank, c] = _gf_mult(M[rank, c], inv, modulus_poly)
            for i in range(rows):
                if i != rank and M[i, j] != 0:
                    f = M[i, j]
                    for c in range(j, cols):
                        M[i, c] = _gf_add(M[i, c], _gf_mult(f, M[rank, c], modulus_poly))
            rank += 1
        return rank

    X_le = np.zeros((le, le), dtype=int)
    Y_me = np.zeros((me, me), dtype=int)
    for i in range(le): X_le[i, (i+1) % le] = 1
    for j in range(me): Y_me[j, (j+1) % me] = 1

    def kron_pow(i, j):
        return np.kron(np.linalg.matrix_power(X_le, i),
                       np.linalg.matrix_power(Y_me, j))

    for alpha in alphas:
        for beta in betas:
            A_sub = np.zeros((le*me, le*me), dtype=object)
            B_sub = np.zeros((le*me, le*me), dtype=object)
            for (i, j) in A_exp:
                coeff = _gf_mult(_power_in_gf(alpha, i, modulus, 1<<r),
                                 _power_in_gf(beta,  j, modulus, 1<<r), modulus)
                T = kron_pow(i, j)
                A_sub[T == 1] = [_gf_add(v, coeff) for v in A_sub[T == 1]]
            for (i, j) in B_exp:
                coeff = _gf_mult(_power_in_gf(alpha, i, modulus, 1<<r),
                                 _power_in_gf(beta,  j, modulus, 1<<r), modulus)
                T = kron_pow(i, j)
                B_sub[T == 1] = [_gf_add(v, coeff) for v in B_sub[T == 1]]
            H_sub = np.concatenate((A_sub, B_sub), axis=1)
            total_rank += rank_gf2r(H_sub, modulus, 1<<r)

    rk_HX = total_rank
    k = n - 2 * rk_HX
    return int(k)

def gf_trace(x, r, mod):
    res, tmp = 0, x
    for _ in range(r):
        res ^= tmp
        tmp = _gf_mult(tmp, tmp, mod)
    return res & 1

def inverse_FT_vec_t(alpha, beta, l, m, r, mod, field_size, t_coeff):
    out = np.empty(l*m, np.uint8)
    t = 0
    for p in range(l):
        a_p = _power_in_gf(alpha, p, mod, field_size)
        for q in range(m):
            b_q = _power_in_gf(beta, q, mod, field_size)
            out[t] = gf_trace(_gf_mult(t_coeff, _gf_mult(a_p, b_q, mod), mod), r, mod)
            t += 1
    return out

def roll_xy(vec, l, m, dx=0, dy=0):
    arr = vec.reshape(l, m)
    if dx: arr = np.roll(arr, dx % l, axis=0)
    if dy: arr = np.roll(arr, dy % m, axis=1)
    return arr.reshape(l*m)

def gf2_rref_basis(vecs, ncols):
    if not vecs:
        return np.zeros((0, ncols), dtype=np.uint8)
    V = (np.array(vecs, dtype=np.uint8) & 1)
    rows, cols = V.shape
    r = c = 0
    while r < rows and c < cols:
        piv = None
        for i in range(r, rows):
            if V[i, c]:
                piv = i; break
        if piv is None:
            c += 1; continue
        if piv != r:
            V[[r, piv]] = V[[piv, r]]
        for i in range(rows):
            if i != r and V[i, c]:
                V[i, :] ^= V[r, :]
        r += 1; c += 1
    return V[:r]

def gf2_rref_with_pivots(M):
    M = (M.copy() & 1).astype(np.uint8)
    rows, cols = M.shape
    r = c = 0
    pivots = []
    while r < rows and c < cols:
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i; break
        if piv is None:
            c += 1; continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        pivots.append(c)
        for i in range(rows):
            if i != r and M[i, c]:
                M[i, :] ^= M[r, :]
        r += 1; c += 1
    return M, pivots

def gf2_nullspace(M):
    M_rref, pivots = gf2_rref_with_pivots(M)
    rows, cols = M_rref.shape
    pivot_set = set(pivots)
    frees = [j for j in range(cols) if j not in pivot_set]
    basis = []
    for f in frees:
        v = np.zeros(cols, dtype=np.uint8)
        v[f] = 1
        for i, p in enumerate(pivots):
            if M_rref[i, f]:
                v[p] ^= 1
        basis.append(v)
    return basis

def extend_basis_to_target(S0, target_basis):
    S = S0.copy()
    if S.size == 0:
        cur_rank = 0
    else:
        cur_rank = gf2_rref_basis([row for row in S], S.shape[1]).shape[0]
    want = len(target_basis)
    if cur_rank >= want:
        return gf2_rref_basis([row for row in S], S.shape[1])
    for v in target_basis:
        cand = np.vstack([S, v[None, :]]) if S.size else v[None, :]
        new_rank = gf2_rref_basis([row for row in cand], cand.shape[1]).shape[0]
        if new_rank > cur_rank:
            S = cand
            cur_rank = new_rank
        if cur_rank >= want:
            break
    return gf2_rref_basis([row for row in S], S.shape[1])

def build_blocks(l, m, A_exp, B_exp):
    Sl = np.eye(l, k=1, dtype=np.uint8); Sl[-1, 0] = 1
    Sm = np.eye(m, k=1, dtype=np.uint8); Sm[-1, 0] = 1
    x = np.kron(Sl, np.eye(m, dtype=np.uint8))
    y = np.kron(np.eye(l, dtype=np.uint8), Sm)
    A_block = sum(np.linalg.matrix_power(x, i) @ np.linalg.matrix_power(y, j) % 2
                  for i, j in A_exp) % 2
    B_block = sum(np.linalg.matrix_power(x, i) @ np.linalg.matrix_power(y, j) % 2
                  for i, j in B_exp) % 2
    HX = np.hstack((A_block, B_block)).astype(np.uint8)
    HZ = np.hstack((B_block.T, A_block.T)).astype(np.uint8)
    Cl = np.zeros((l, l), dtype=np.uint8)
    for j in range(l):
        Cl[(-j) % l, j] = 1
    Cm = np.zeros((m, m), dtype=np.uint8)
    for j in range(m):
        Cm[(-j) % m, j] = 1
    C = np.kron(Cl, Cm).astype(np.uint8)
    return A_block, B_block, HX, HZ, C

def kernel_basis_ifft(l, m, A_exp, B_exp, do_half_period_shifts=True):
    l1, _ = two_pow_factor(l)
    m1, _ = two_pow_factor(m)
    r = max(1, min_ext_deg(l, m))
    mod = (0, 0) if r == 1 else _find_irreducible_polynomial(r)
    field_size = 1 << r
    elems  = [1] if r == 1 else list(range(1, field_size))
    alphas = [1] if l1 == 1 else [a for a in elems if _power_in_gf(a, l1, mod, field_size) == 1]
    betas  = [1] if m1 == 1 else [b for b in elems if _power_in_gf(b, m1, mod, field_size) == 1]
    dx = (l & -l) >> 1
    dy = (m & -m) >> 1
    t_basis = [1 << k for k in range(r)]
    vecs = []
    for a in alphas:
        for b in betas:
            fA = fB = 0
            for i, j in A_exp:
                fA = _gf_add(fA, _gf_mult(_power_in_gf(a, i, mod, field_size),
                                          _power_in_gf(b, j, mod, field_size), mod))
            for i, j in B_exp:
                fB = _gf_add(fB, _gf_mult(_power_in_gf(a, i, mod, field_size),
                                          _power_in_gf(b, j, mod, field_size), mod))
            if fA or fB:
                continue
            for t_coeff in t_basis:
                base = inverse_FT_vec_t(a, b, l, m, r, mod, field_size, t_coeff)
                vecs.append(base)
                if do_half_period_shifts:
                    if dx: vecs.append(roll_xy(base, l, m, dx=dx, dy=0))
                    if dy: vecs.append(roll_xy(base, l, m, dx=0,  dy=dy))
                    if dx and dy: vecs.append(roll_xy(base, l, m, dx=dx, dy=dy))
    S0 = gf2_rref_basis(vecs, l*m)
    A_block, B_block, _, _, _ = build_blocks(l, m, A_exp, B_exp)
    M = np.vstack([A_block, B_block]) % 2
    target_basis = gf2_nullspace(M)
    S = extend_basis_to_target(S0, target_basis)
    return S

def logical_ops_ifft(l, m, A_exp, B_exp):
    S = kernel_basis_ifft(l, m, A_exp, B_exp, do_half_period_shifts=True)
    _, _, HX, HZ, C = build_blocks(l, m, A_exp, B_exp)
    Z_ops = np.hstack([S, np.zeros_like(S, dtype=np.uint8)])
    X_ops = np.hstack([np.zeros_like(S, dtype=np.uint8), (S @ C) % 2])
    return Z_ops.astype(np.uint8), X_ops.astype(np.uint8), HX, HZ

if __name__ == "__main__":
    l_val, m_val = 21, 6*3
    A_exponents = [(3,0),(0,10),(0,17)]
    B_exponents = [(0,5),(3,0),(19,0)]
    k_val = k_bb_code_final(l_val, m_val, A_exponents, B_exponents)  # compute k
    print("k =", k_val)
    Zs, Xs, HX, HZ = logical_ops_ifft(l_val, m_val, A_exponents, B_exponents)  # logical ops
    print("Z_ops", Zs.shape, "X_ops", Xs.shape)
    z_ok = not (HX @ Zs.T % 2).any()
    x_ok = not (HZ @ Xs.T % 2).any()
    print("HX·Z=0 ?", z_ok, "   HZ·X=0 ?", x_ok)
    split = l_val * m_val
    for i in range(Zs.shape[0]):
        print(f"\nZ{i+1}: {''.join(map(str, Zs[i][:split]))} | {''.join(map(str, Zs[i][split:]))}")
        print(f"X{i+1}: {''.join(map(str, Xs[i][:split]))} | {''.join(map(str, Xs[i][split:]))}")
    assert Zs.shape[0] == k_val // 2 == Xs.shape[0]
    assert z_ok and x_ok
