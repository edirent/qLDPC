"""
Fast Fourier-based utilities for bivariate-bicycle qLDPC code calculations.

Main changes relative to the original pasted script:
  * negative and oversized torus exponents are normalized modulo (l, m),
  * duplicate monomials cancel over GF(2),
  * the GF(2^7) modulus bug is removed and irreducible moduli can be generated,
  * dense shift blocks are built directly instead of using large matrix powers,
  * Fourier-rank computation reuses local shift masks and integer arrays,
  * a canonical CSS logical basis constructor is provided,
  * Frobenius-orbit compression and specialized finite-field kernels speed up k computation.

Elements of GF(2^r) are represented as integers whose bits are polynomial
coefficients over GF(2). Addition is XOR; multiplication reduces modulo a
monic irreducible polynomial.
"""

from __future__ import annotations

from functools import lru_cache
from math import gcd
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

Exponent = Tuple[int, int]

# ---------------------------------------------------------------------------
# Exponent and dimension handling
# ---------------------------------------------------------------------------


def _validate_lattice(l: int, m: int) -> Tuple[int, int]:
    """Validate and return positive integer lattice dimensions."""
    l = int(l)
    m = int(m)
    if l <= 0 or m <= 0:
        raise ValueError(f"l and m must be positive; got l={l}, m={m}.")
    return l, m


def normalize_exponents(
    exps: Iterable[Sequence[int]],
    l: int,
    m: int,
    *,
    cancel_duplicates: bool = True,
) -> List[Exponent]:
    """Normalize monomial exponents on the l-by-m torus.

    BB-code polynomials are over F2[x, y]/<x^l-1, y^m-1>. Therefore negative
    exponents and exponents outside the fundamental domain are reduced modulo
    l and m. If ``cancel_duplicates`` is true, repeated monomials cancel as
    they should over GF(2).
    """
    l, m = _validate_lattice(l, m)
    if cancel_duplicates:
        parity: Dict[Exponent, int] = {}
        for pair in exps:
            if len(pair) != 2:
                raise ValueError(f"Expected exponent pair (i, j); got {pair!r}.")
            key = (int(pair[0]) % l, int(pair[1]) % m)
            parity[key] = parity.get(key, 0) ^ 1
        return sorted([key for key, bit in parity.items() if bit])

    out: List[Exponent] = []
    for pair in exps:
        if len(pair) != 2:
            raise ValueError(f"Expected exponent pair (i, j); got {pair!r}.")
        out.append((int(pair[0]) % l, int(pair[1]) % m))
    return out


# ---------------------------------------------------------------------------
# Polynomial arithmetic over GF(2), used to find irreducible moduli
# ---------------------------------------------------------------------------


def _poly_degree(poly: int) -> int:
    """Return the highest nonzero coefficient index of an encoded polynomial."""
    return int(poly).bit_length() - 1


def _poly_mod(poly: int, modulus: int) -> int:
    """Return poly modulo modulus for binary polynomials encoded as integers."""
    if modulus <= 0:
        raise ValueError("Polynomial modulus must be positive.")
    mod_degree = _poly_degree(modulus)
    poly = int(poly)
    while poly and _poly_degree(poly) >= mod_degree:
        poly ^= modulus << (_poly_degree(poly) - mod_degree)
    return poly


def _poly_gcd(a: int, b: int) -> int:
    """GCD of two GF(2) polynomials encoded as integers."""
    a, b = int(a), int(b)
    while b:
        a, b = b, _poly_mod(a, b)
    return a


def _poly_mul_mod(a: int, b: int, modulus: int) -> int:
    """Multiply two GF(2) polynomials and reduce modulo modulus."""
    res = 0
    a = int(a)
    b = int(b)
    while b:
        if b & 1:
            res ^= a
        a <<= 1
        b >>= 1
    return _poly_mod(res, modulus)


def _poly_square_mod(a: int, modulus: int) -> int:
    """Square a GF(2) polynomial and reduce it modulo ``modulus``."""
    return _poly_mul_mod(a, a, modulus)


def _is_irreducible_poly(poly: int, degree: int) -> bool:
    """Rabin-style irreducibility test over GF(2)."""
    if degree <= 0:
        return False
    if _poly_degree(poly) != degree:
        return False
    if not (poly & 1):
        return False
    if degree == 1:
        return poly == 0b11

    x = 0b10

    # Check x^(2^degree) == x mod poly.
    xp = x
    for _ in range(degree):
        xp = _poly_square_mod(xp, poly)
    if xp != x:
        return False

    # Check gcd(x^(2^i)-x, poly)=1 for i=1..floor(degree/2).
    xp = x
    for _ in range(1, degree // 2 + 1):
        xp = _poly_square_mod(xp, poly)
        if _poly_gcd(xp ^ x, poly) != 1:
            return False
    return True


@lru_cache(maxsize=None)
def _find_irreducible_polynomial(degree: int) -> Tuple[int, int]:
    """Return a monic irreducible polynomial of the requested degree.

    The returned pair is ``(poly, high_bit)``, where ``high_bit == 1 << degree``.
    The previous script used a reducible degree-7 polynomial; this implementation
    verifies the table entries and falls back to a deterministic search.
    """
    degree = int(degree)
    if degree < 1:
        raise ValueError(f"Extension degree must be at least 1; got {degree}.")

    preferred = {
        1: 0b11,          # x + 1
        2: 0b111,         # x^2 + x + 1
        3: 0b1011,        # x^3 + x + 1
        4: 0b10011,       # x^4 + x + 1
        5: 0b100101,      # x^5 + x^2 + 1
        6: 0b1000011,     # x^6 + x + 1
        7: 0b10000011,    # x^7 + x + 1
        8: 0b100011101,   # x^8 + x^4 + x^3 + x^2 + 1
        9: 0b1000010001,  # x^9 + x^4 + 1
        10: 0b10000001001,
        11: 0b100000000101,
        12: 0b1000001010011,
        13: 0b10000000011011,
        14: 0b100000000101011,
        15: 0b1000000000000011,
        16: 0b10001000000001011,
    }
    candidate = preferred.get(degree)
    if candidate is not None and _is_irreducible_poly(candidate, degree):
        return candidate, 1 << degree

    # Deterministic search over monic polynomials with nonzero constant term.
    start = (1 << degree) | 1
    stop = 1 << (degree + 1)
    for poly in range(start, stop, 2):
        if _is_irreducible_poly(poly, degree):
            return poly, 1 << degree
    raise ValueError(f"Could not find an irreducible polynomial of degree {degree}.")


# ---------------------------------------------------------------------------
# Basic arithmetic over GF(2^r)
# ---------------------------------------------------------------------------


def _gf_mult(a: int, b: int, modulus: Tuple[int, int]) -> int:
    """Multiply two GF(2^r) elements encoded as integers."""
    a = int(a)
    b = int(b)
    poly, high_bit = modulus
    if high_bit == 0:  # GF(2), with elements restricted to 0 or 1.
        return (a & 1) & (b & 1)

    res = 0
    while b > 0:
        if b & 1:
            res ^= a
        b >>= 1
        a <<= 1
        if a & high_bit:
            a ^= poly
    return res


def _gf_add(a: int, b: int) -> int:
    """Add two field elements; in characteristic two this is bitwise XOR."""
    return int(a) ^ int(b)


def _gf_inverse(a: int, modulus: Tuple[int, int], field_size: int) -> int:
    """Return a^{-1} in GF(2^r) by Fermat exponentiation."""
    a = int(a)
    if a == 0:
        raise ValueError("Cannot compute inverse of 0.")
    return _power_in_gf(a, field_size - 2, modulus, field_size)


@lru_cache(maxsize=None)
def _gf_mult_table(modulus_poly: int, high_bit: int, field_size: int) -> np.ndarray:
    """Multiplication table for GF(2^r), used by vectorized elimination."""
    modulus = (int(modulus_poly), int(high_bit))
    table = np.zeros((int(field_size), int(field_size)), dtype=np.int64)
    for a in range(int(field_size)):
        for b in range(int(field_size)):
            table[a, b] = _gf_mult(a, b, modulus)
    return table


def _power_in_gf(base: int, exp: int, modulus: Tuple[int, int], field_size: int) -> int:
    """Exponentiate a field element by repeated squaring.

    Negative exponents are supported for nonzero bases.
    """
    base = int(base)
    exp = int(exp)
    if exp < 0:
        if base == 0:
            raise ValueError("Cannot raise 0 to a negative power in a field.")
        base = _gf_inverse(base, modulus, field_size)
        exp = -exp
    if exp == 0:
        return 1

    res = 1
    cur = base
    while exp > 0:
        if exp & 1:
            res = _gf_mult(res, cur, modulus)
        cur = _gf_mult(cur, cur, modulus)
        exp >>= 1
    return res


# ---------------------------------------------------------------------------
# Lattice factorization and roots of unity
# ---------------------------------------------------------------------------


def _get_ord_2(n: int) -> int:
    """Return the multiplicative order of 2 modulo odd n."""
    n = int(n)
    if n == 1:
        return 1
    if n <= 0 or gcd(2, n) != 1:
        raise ValueError(f"n must be a positive odd integer; got {n}.")
    k, res = 1, 2 % n
    while res != 1:
        res = (res * 2) % n
        k += 1
    return k


def two_pow_factor(n: int) -> Tuple[int, int]:
    """Split n as n = odd * twop, where twop is a power of two."""
    n = int(n)
    if n <= 0:
        raise ValueError(f"n must be positive; got {n}.")
    odd, twop = n, 1
    while odd % 2 == 0:
        odd //= 2
        twop *= 2
    return odd, twop


def min_ext_deg(l: int, m: int) -> int:
    """Minimum extension degree containing roots for the odd parts of l and m."""
    l, m = _validate_lattice(l, m)
    l1, _ = two_pow_factor(l)
    m1, _ = two_pow_factor(m)
    if l1 == m1 == 1:
        return 1
    return int(np.lcm(_get_ord_2(l1), _get_ord_2(m1)))


def _field_context(l: int, m: int):
    """Return factorization, finite-field context, and root lists for (l, m)."""
    l, m = _validate_lattice(l, m)
    l1, le = two_pow_factor(l)
    m1, me = two_pow_factor(m)
    r = min_ext_deg(l, m)

    if r == 1:
        modulus = (0, 0)
        field_size = 2
        elems = [1]
    else:
        modulus = _find_irreducible_polynomial(r)
        field_size = 1 << r
        elems = list(range(1, field_size))

    def roots(order: int) -> List[int]:
        """Enumerate all ``order``-th roots of unity in the selected field."""
        if order == 1:
            return [1]
        out = [a for a in elems if _power_in_gf(a, order, modulus, field_size) == 1]
        if len(out) != order:
            raise RuntimeError(
                f"Expected {order} roots of unity in GF(2^{r}); found {len(out)}."
            )
        return out

    return l1, le, m1, me, r, modulus, field_size, roots(l1), roots(m1)


# ---------------------------------------------------------------------------
# GF(2) linear algebra
# ---------------------------------------------------------------------------


class _GF2RowSpace:
    """Incremental row-space membership over GF(2)."""

    def __init__(self, ncols: int):
        """Create an empty row space with rows of fixed width ``ncols``."""
        self.ncols = int(ncols)
        self._rows_by_pivot: Dict[int, np.ndarray] = {}

    def add(self, row: np.ndarray) -> bool:
        """Add row if independent. Return True iff rank increased."""
        v = (np.asarray(row, dtype=np.uint8).copy() & 1)
        if v.size != self.ncols:
            raise ValueError(f"Expected row of length {self.ncols}; got {v.size}.")
        while True:
            nz = np.flatnonzero(v)
            if nz.size == 0:
                return False
            pivot = int(nz[0])
            existing = self._rows_by_pivot.get(pivot)
            if existing is None:
                self._rows_by_pivot[pivot] = v
                return True
            v ^= existing

    @property
    def rank(self) -> int:
        """Current number of independent rows in the incremental basis."""
        return len(self._rows_by_pivot)

    def basis_matrix(self) -> np.ndarray:
        """Return the stored pivot rows as a binary matrix sorted by pivot."""
        if not self._rows_by_pivot:
            return np.zeros((0, self.ncols), dtype=np.uint8)
        return np.vstack([self._rows_by_pivot[p] for p in sorted(self._rows_by_pivot)])


def _as_binary_matrix(M: Iterable, ncols: int | None = None) -> np.ndarray:
    """Convert row-like input to a 2D uint8 matrix with entries reduced mod 2.

    A one-dimensional input is treated as one row.  Empty non-array iterables
    need ``ncols`` so that downstream linear-algebra code still receives a
    correctly shaped ``(0, ncols)`` matrix.
    """
    if isinstance(M, np.ndarray):
        A = (M.copy() & 1).astype(np.uint8)
        if A.ndim == 1:
            if ncols is None:
                ncols = A.size
            A = A.reshape(1, ncols)
        return A
    rows = list(M)
    if not rows:
        if ncols is None:
            raise ValueError("ncols is required for an empty matrix/list.")
        return np.zeros((0, ncols), dtype=np.uint8)
    return (np.array(rows, dtype=np.uint8) & 1)


def gf2_rref_basis(vecs: Iterable, ncols: int | None = None) -> np.ndarray:
    """Return a reduced independent row basis for binary vectors."""
    V = _as_binary_matrix(vecs, ncols)
    if V.size == 0:
        return V
    rows, cols = V.shape
    r = c = 0
    while r < rows and c < cols:
        piv = None
        for i in range(r, rows):
            if V[i, c]:
                piv = i
                break
        if piv is None:
            c += 1
            continue
        if piv != r:
            V[[r, piv]] = V[[piv, r]]
        for i in range(rows):
            if i != r and V[i, c]:
                V[i, :] ^= V[r, :]
        r += 1
        c += 1
    return V[:r]


def gf2_rref_with_pivots(M: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    """Return binary RREF of M together with pivot column indices."""
    M = (np.asarray(M, dtype=np.uint8).copy() & 1)
    rows, cols = M.shape
    r = c = 0
    pivots: List[int] = []
    while r < rows and c < cols:
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            c += 1
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        pivots.append(c)
        for i in range(rows):
            if i != r and M[i, c]:
                M[i, :] ^= M[r, :]
        r += 1
        c += 1
    return M, pivots


def gf2_rank(M: np.ndarray) -> int:
    """Rank of a binary matrix."""
    _, pivots = gf2_rref_with_pivots(M)
    return len(pivots)


def gf2_nullspace(M: np.ndarray) -> np.ndarray:
    """Compute a binary nullspace basis as rows."""
    M_rref, pivots = gf2_rref_with_pivots(M)
    _, cols = M_rref.shape
    pivot_set = set(pivots)
    frees = [j for j in range(cols) if j not in pivot_set]
    if not frees:
        return np.zeros((0, cols), dtype=np.uint8)

    basis = []
    for f in frees:
        v = np.zeros(cols, dtype=np.uint8)
        v[f] = 1
        for i, p in enumerate(pivots):
            if M_rref[i, f]:
                v[p] ^= 1
        basis.append(v)
    return np.vstack(basis).astype(np.uint8)


def gf2_inverse_square(M: np.ndarray) -> np.ndarray:
    """Inverse of a nonsingular square binary matrix."""
    M = (np.asarray(M, dtype=np.uint8).copy() & 1)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("M must be square.")
    n = M.shape[0]
    aug = np.hstack([M, np.eye(n, dtype=np.uint8)])
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, n):
            if aug[i, c]:
                piv = i
                break
        if piv is None:
            raise ValueError("Matrix is singular over GF(2).")
        if piv != r:
            aug[[r, piv]] = aug[[piv, r]]
        for i in range(n):
            if i != r and aug[i, c]:
                aug[i, :] ^= aug[r, :]
        r += 1
    return aug[:, n:].astype(np.uint8)


def gf2_matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Matrix product over GF(2)."""
    return ((np.asarray(A, dtype=np.uint8) @ np.asarray(B, dtype=np.uint8)) & 1).astype(np.uint8)


def _quotient_basis(space_basis: np.ndarray, mod_basis: np.ndarray, ncols: int) -> np.ndarray:
    """Return a basis of span(space_basis)/span(mod_basis), represented by rows."""
    rowspace = _GF2RowSpace(ncols)
    for row in _as_binary_matrix(mod_basis, ncols):
        rowspace.add(row)

    chosen: List[np.ndarray] = []
    for row in _as_binary_matrix(space_basis, ncols):
        before = rowspace.rank
        if rowspace.add(row) and rowspace.rank > before:
            chosen.append((row.copy() & 1).astype(np.uint8))
    if not chosen:
        return np.zeros((0, ncols), dtype=np.uint8)
    return np.vstack(chosen).astype(np.uint8)


def extend_basis_to_target(S0: np.ndarray, target_basis: np.ndarray, ncols: int | None = None) -> np.ndarray:
    """Extend an initial independent set until it spans a target space."""
    if ncols is None:
        if S0.size:
            ncols = S0.shape[1]
        elif target_basis.size:
            ncols = target_basis.shape[1]
        else:
            raise ValueError("ncols is needed when both inputs are empty.")
    ncols = int(ncols)

    rowspace = _GF2RowSpace(ncols)
    chosen: List[np.ndarray] = []
    for row in _as_binary_matrix(S0, ncols):
        if rowspace.add(row):
            chosen.append(row.copy())

    want = _as_binary_matrix(target_basis, ncols).shape[0]
    for row in _as_binary_matrix(target_basis, ncols):
        if len(chosen) >= want:
            break
        if rowspace.add(row):
            chosen.append(row.copy())

    if not chosen:
        return np.zeros((0, ncols), dtype=np.uint8)
    return gf2_rref_basis(np.vstack(chosen), ncols)


# ---------------------------------------------------------------------------
# Dense shift blocks and BB code blocks
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _shift_matrix_cached(l: int, m: int, dx: int, dy: int) -> np.ndarray:
    """Permutation matrix for x^dx y^dy on a flattened l-by-m grid."""
    l, m = _validate_lattice(l, m)
    dx %= l
    dy %= m
    n0 = l * m
    rows = np.arange(n0, dtype=np.int64)
    p = rows // m
    q = rows % m
    cols = ((p + dx) % l) * m + ((q + dy) % m)
    out = np.zeros((n0, n0), dtype=np.uint8)
    out[rows, cols] = 1
    return out


def _sum_shift_blocks(l: int, m: int, exps: Sequence[Exponent]) -> np.ndarray:
    """XOR the torus shift matrices for all monomials in a BB polynomial."""
    n0 = l * m
    block = np.zeros((n0, n0), dtype=np.uint8)
    for i, j in exps:
        # Exponents are already parity-normalized; XOR remains correct.
        block ^= _shift_matrix_cached(l, m, i, j)
    return block


def _coordinate_inversion_matrix(l: int, m: int) -> np.ndarray:
    """Permutation matrix mapping each torus coordinate ``(p, q)`` to ``(-p, -q)``."""
    n0 = l * m
    cols = np.arange(n0, dtype=np.int64)
    p = cols // m
    q = cols % m
    rows = ((-p) % l) * m + ((-q) % m)
    C = np.zeros((n0, n0), dtype=np.uint8)
    C[rows, cols] = 1
    return C


def build_blocks(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]):
    """Construct A, B, H_X, H_Z, and the coordinate-inversion matrix C."""
    l, m = _validate_lattice(l, m)
    A_exp = normalize_exponents(A_exp, l, m)
    B_exp = normalize_exponents(B_exp, l, m)

    A_block = _sum_shift_blocks(l, m, A_exp)
    B_block = _sum_shift_blocks(l, m, B_exp)
    HX = np.hstack((A_block, B_block)).astype(np.uint8)
    HZ = np.hstack((B_block.T, A_block.T)).astype(np.uint8)
    C = _coordinate_inversion_matrix(l, m)
    return A_block, B_block, HX, HZ, C


# ---------------------------------------------------------------------------
# Rank computation by Fourier block decomposition
# ---------------------------------------------------------------------------


def _rank_gf2r(matrix: np.ndarray, modulus: Tuple[int, int], field_size: int) -> int:
    """Rank of a matrix over GF(2^r).

    Elimination is vectorized with a precomputed multiplication table. Addition
    in characteristic two is XOR, so eliminating a row is a table lookup followed
    by one vectorized XOR.
    """
    M = np.asarray(matrix, dtype=np.int64).copy()
    rows, cols = M.shape
    rank = 0
    mult_table = _gf_mult_table(int(modulus[0]), int(modulus[1]), int(field_size))

    for j in range(cols):
        nonzero = np.flatnonzero(M[rank:, j])
        if nonzero.size == 0:
            continue
        piv = rank + int(nonzero[0])
        if piv != rank:
            M[[rank, piv]] = M[[piv, rank]]

        pv = int(M[rank, j])
        inv = _gf_inverse(pv, modulus, field_size)
        if inv != 1:
            M[rank, j:] = mult_table[inv, M[rank, j:]]

        elim_rows = np.flatnonzero(M[:, j])
        for i in elim_rows:
            i = int(i)
            if i == rank:
                continue
            factor = int(M[i, j])
            M[i, j:] ^= mult_table[factor, M[rank, j:]]

        rank += 1
        if rank == rows:
            break
    return rank


# ---------------------------------------------------------------------------
# Optimized Fourier-rank core
# ---------------------------------------------------------------------------


def _prime_factors_unique(n: int) -> Tuple[int, ...]:
    """Return the distinct prime factors of n."""
    n = int(n)
    factors: List[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors.append(n)
    return tuple(factors)


@lru_cache(maxsize=None)
def _gf_tables_fast(modulus_poly: int, high_bit: int, field_size: int):
    """Finite-field multiplication/inversion tables using a primitive element."""
    modulus_poly = int(modulus_poly)
    high_bit = int(high_bit)
    field_size = int(field_size)
    if field_size == 2:
        return ((0, 0), (0, 1)), (0, 1), (1,), (0, 0)

    modulus = (modulus_poly, high_bit)
    order = field_size - 1
    prime_factors = _prime_factors_unique(order)
    primitive = None
    for candidate in range(2, field_size):
        if all(_power_in_gf(candidate, order // p, modulus, field_size) != 1 for p in prime_factors):
            primitive = candidate
            break
    if primitive is None:
        raise RuntimeError(f"Could not find primitive element for GF({field_size}).")

    exp_arr = np.empty(order, dtype=np.int64)
    log_arr = np.zeros(field_size, dtype=np.int64)
    value = 1
    for exponent in range(order):
        exp_arr[exponent] = value
        log_arr[value] = exponent
        value = _gf_mult(value, primitive, modulus)
    if value != 1:
        raise RuntimeError("Primitive-element table construction failed.")

    idx = (log_arr[:, None] + log_arr[None, :]) % order
    mult_arr = exp_arr[idx]
    mult_arr[0, :] = 0
    mult_arr[:, 0] = 0

    inv_arr = np.zeros(field_size, dtype=np.int64)
    inv_arr[1:] = exp_arr[(order - log_arr[1:]) % order]
    return (
        tuple(tuple(int(x) for x in row) for row in mult_arr.tolist()),
        tuple(int(x) for x in inv_arr.tolist()),
        tuple(int(x) for x in exp_arr.tolist()),
        tuple(int(x) for x in log_arr.tolist()),
    )


@lru_cache(maxsize=None)
def _field_context_fast(l: int, m: int):
    """Fast field context and root lists for Fourier decomposition."""
    l, m = _validate_lattice(l, m)
    l1, le = two_pow_factor(l)
    m1, me = two_pow_factor(m)
    r = min_ext_deg(l, m)
    if r == 1:
        modulus = (0, 0)
        field_size = 2
    else:
        modulus = _find_irreducible_polynomial(r)
        field_size = 1 << r
    mult_table, inverse_table, exp_table, log_table = _gf_tables_fast(modulus[0], modulus[1], field_size)

    def roots(order: int) -> Tuple[int, ...]:
        """Return roots of unity using the cached primitive-element tables."""
        order = int(order)
        if order == 1:
            return (1,)
        group_order = field_size - 1
        if group_order % order != 0:
            raise RuntimeError(f"Expected root order {order} to divide {group_order}.")
        step = group_order // order
        return tuple(int(exp_table[(j * step) % group_order]) for j in range(order))

    alphas = roots(l1)
    betas = roots(m1)
    return l1, le, m1, me, r, modulus, field_size, alphas, betas, mult_table, inverse_table, exp_table, log_table


@lru_cache(maxsize=None)
def _local_shift_permutations(le: int, me: int) -> Dict[Exponent, Tuple[int, ...]]:
    """Local power-of-two shift permutations for repeated-root blocks."""
    le = int(le)
    me = int(me)
    if le <= 0 or me <= 0:
        raise ValueError("le and me must be positive.")
    block_dim = le * me
    out: Dict[Exponent, Tuple[int, ...]] = {}
    for dx in range(le):
        for dy in range(me):
            cols: List[int] = []
            for row in range(block_dim):
                p = row // me
                q = row % me
                cols.append(((p + dx) % le) * me + ((q + dy) % me))
            out[(dx, dy)] = tuple(cols)
    return out


@lru_cache(maxsize=None)
def _frobenius_pair_orbits(l: int, m: int) -> Tuple[Tuple[int, int, int], ...]:
    """Representatives and sizes of Frobenius orbits of root pairs."""
    _, _, _, _, _, _, _, alphas, betas, mult_table, _, _, _ = _field_context_fast(l, m)
    unseen = {(a, b) for a in alphas for b in betas}
    reps: List[Tuple[int, int, int]] = []
    while unseen:
        alpha, beta = next(iter(unseen))
        orbit: List[Tuple[int, int]] = []
        a, b = alpha, beta
        while (a, b) not in orbit:
            orbit.append((a, b))
            a = mult_table[a][a]
            b = mult_table[b][b]
        for pair in orbit:
            unseen.discard(pair)
        reps.append((alpha, beta, len(orbit)))
    return tuple(reps)


def _field_power_from_log(a: int, exponent: int, field_size: int, exp_table: Tuple[int, ...], log_table: Tuple[int, ...]) -> int:
    """Power of a nonzero field element using log/antilog tables."""
    if field_size == 2:
        return 1
    order = field_size - 1
    return exp_table[(log_table[int(a)] * (int(exponent) % order)) % order]


def _rank_gf2r_rows_inplace(rows: List[List[int]], mult_table, inverse_table) -> int:
    """Rank over GF(2^r) for a tiny dense row matrix, in place."""
    nrows = len(rows)
    if nrows == 0:
        return 0
    ncols = len(rows[0])
    rank = 0
    for col in range(ncols):
        pivot = -1
        for rr in range(rank, nrows):
            if rows[rr][col]:
                pivot = rr
                break
        if pivot < 0:
            continue
        if pivot != rank:
            rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = inverse_table[rows[rank][col]]
        if inv != 1:
            scale = mult_table[inv]
            prow = rows[rank]
            for cc in range(col, ncols):
                prow[cc] = scale[prow[cc]]
        prow = rows[rank]
        for rr in range(nrows):
            if rr != rank and rows[rr][col]:
                scale = mult_table[rows[rr][col]]
                row = rows[rr]
                for cc in range(col, ncols):
                    row[cc] ^= scale[prow[cc]]
        rank += 1
        if rank == nrows:
            break
    return rank


def _rank_gf2r_two_rows(rows: List[List[int]], mult_table, inverse_table) -> int:
    """Specialized rank for a 2 x c matrix over GF(2^r)."""
    r0 = rows[0]
    r1 = rows[1]
    first0 = next((i for i, v in enumerate(r0) if v), -1)
    first1 = next((i for i, v in enumerate(r1) if v), -1)
    if first0 < 0:
        return 1 if first1 >= 0 else 0
    if first1 < 0:
        return 1
    lam = mult_table[r1[first0]][inverse_table[r0[first0]]] if r1[first0] else 0
    scale = mult_table[lam]
    for idx, val in enumerate(r0):
        if r1[idx] != scale[val]:
            return 2
    return 1


def k_bb_code_final(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> int:
    """Compute k for a BB-style CSS code using optimized Fourier decomposition.

    Optimizations:
      * root pairs are quotiented by Frobenius orbits;
      * finite-field multiplication, inverses, and powers use log/antilog tables;
      * local repeated-root blocks are assembled row-wise from permutations;
      * 1-row and 2-row rank computations use specialized paths.
    """
    l, m = _validate_lattice(l, m)
    A_exp = tuple(normalize_exponents(A_exp, l, m))
    B_exp = tuple(normalize_exponents(B_exp, l, m))
    n = 2 * l * m

    _, le, _, me, _, _, field_size, _, _, mult_table, inverse_table, exp_table, log_table = _field_context_fast(l, m)
    block_dim = le * me
    local_perms = _local_shift_permutations(le, me)
    exps = tuple(sorted(set(A_exp + B_exp)))

    total_rank = 0
    for alpha, beta, orbit_size in _frobenius_pair_orbits(l, m):
        alpha_pows = {i: _field_power_from_log(alpha, i, field_size, exp_table, log_table) for i, _ in exps}
        beta_pows = {j: _field_power_from_log(beta, j, field_size, exp_table, log_table) for _, j in exps}

        # Residue-field/unit shortcut.  In the repeated-root local ring, a local
        # polynomial is a unit iff its residue value is nonzero.  If A or B is a
        # unit, the row block [A B] is surjective and has full row rank.  Only
        # common-zero blocks require local Gaussian elimination.
        a_val = 0
        b_val = 0
        for i, j in A_exp:
            a_val ^= mult_table[alpha_pows[i]][beta_pows[j]]
        for i, j in B_exp:
            b_val ^= mult_table[alpha_pows[i]][beta_pows[j]]
        if a_val or b_val:
            total_rank += orbit_size * block_dim
            continue

        if block_dim == 1:
            # Both residue values vanish and there is no repeated-root local
            # part, so the entire 1x2 block is zero.
            continue

        rows = [[0] * (2 * block_dim) for _ in range(block_dim)]
        for i, j in A_exp:
            coeff = mult_table[alpha_pows[i]][beta_pows[j]]
            if coeff:
                perm = local_perms[(i % le, j % me)]
                for row, col in enumerate(perm):
                    rows[row][col] ^= coeff
        for i, j in B_exp:
            coeff = mult_table[alpha_pows[i]][beta_pows[j]]
            if coeff:
                perm = local_perms[(i % le, j % me)]
                offset = block_dim
                for row, col in enumerate(perm):
                    rows[row][offset + col] ^= coeff

        if block_dim == 2:
            rank = _rank_gf2r_two_rows(rows, mult_table, inverse_table)
        else:
            rank = _rank_gf2r_rows_inplace(rows, mult_table, inverse_table)
        total_rank += orbit_size * rank

    return int(n - 2 * total_rank)

def direct_k_bb_code(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> int:
    """Compute k directly from GF(2) ranks of H_X and H_Z."""
    _, _, HX, HZ, _ = build_blocks(l, m, A_exp, B_exp)
    n = HX.shape[1]
    return int(n - gf2_rank(HX) - gf2_rank(HZ))


# ---------------------------------------------------------------------------
# Inverse Fourier reconstruction of ker(A) intersect ker(B)
# ---------------------------------------------------------------------------


def gf_trace(x: int, r: int, mod: Tuple[int, int]) -> int:
    """Trace GF(2^r) -> GF(2): Tr(x)=x+x^2+...+x^(2^(r-1))."""
    if r == 1:
        return int(x) & 1
    res, tmp = 0, int(x)
    for _ in range(r):
        res ^= tmp
        tmp = _gf_mult(tmp, tmp, mod)
    return res & 1


def inverse_FT_vec_t(
    alpha: int,
    beta: int,
    l: int,
    m: int,
    r: int,
    mod: Tuple[int, int],
    field_size: int,
    t_coeff: int,
) -> np.ndarray:
    """Build one binary trace vector from a Fourier mode and trace coefficient."""
    alpha_pows = [_power_in_gf(alpha, p, mod, field_size) for p in range(l)]
    beta_pows = [_power_in_gf(beta, q, mod, field_size) for q in range(m)]
    out = np.empty(l * m, dtype=np.uint8)
    t = 0
    for p in range(l):
        a_p = alpha_pows[p]
        for q in range(m):
            out[t] = gf_trace(_gf_mult(t_coeff, _gf_mult(a_p, beta_pows[q], mod), mod), r, mod)
            t += 1
    return out


def roll_xy(vec: np.ndarray, l: int, m: int, dx: int = 0, dy: int = 0) -> np.ndarray:
    """Cyclically shift a flattened l-by-m binary vector by (dx, dy)."""
    arr = np.asarray(vec, dtype=np.uint8).reshape(l, m)
    if dx:
        arr = np.roll(arr, dx % l, axis=0)
    if dy:
        arr = np.roll(arr, dy % m, axis=1)
    return arr.reshape(l * m).copy()


def kernel_basis_ifft(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    do_half_period_shifts: bool = True,
) -> np.ndarray:
    """Recover a binary basis for ker(A) intersect ker(B).

    Inverse-Fourier candidates are generated first and then completed against
    the direct GF(2) nullspace of [A; B], so the returned basis is exact even
    when the candidate set is incomplete for repeated-root/even-period cases.
    """
    l, m = _validate_lattice(l, m)
    A_exp = normalize_exponents(A_exp, l, m)
    B_exp = normalize_exponents(B_exp, l, m)
    l1, _, m1, _, r, mod, field_size, alphas, betas = _field_context(l, m)

    dx = (l & -l) >> 1
    dy = (m & -m) >> 1
    t_basis = [1 << k for k in range(r)]

    vecs: List[np.ndarray] = []
    for a in alphas:
        for b in betas:
            fA = 0
            fB = 0
            for i, j in A_exp:
                fA = _gf_add(
                    fA,
                    _gf_mult(_power_in_gf(a, i, mod, field_size), _power_in_gf(b, j, mod, field_size), mod),
                )
            for i, j in B_exp:
                fB = _gf_add(
                    fB,
                    _gf_mult(_power_in_gf(a, i, mod, field_size), _power_in_gf(b, j, mod, field_size), mod),
                )
            if fA or fB:
                continue
            for t_coeff in t_basis:
                base = inverse_FT_vec_t(a, b, l, m, r, mod, field_size, t_coeff)
                vecs.append(base)
                if do_half_period_shifts:
                    if dx:
                        vecs.append(roll_xy(base, l, m, dx=dx, dy=0))
                    if dy:
                        vecs.append(roll_xy(base, l, m, dx=0, dy=dy))
                    if dx and dy:
                        vecs.append(roll_xy(base, l, m, dx=dx, dy=dy))

    S0 = gf2_rref_basis(vecs, l * m) if vecs else np.zeros((0, l * m), dtype=np.uint8)
    A_block, B_block, _, _, _ = build_blocks(l, m, A_exp, B_exp)
    target_basis = gf2_nullspace(np.vstack([A_block, B_block]) % 2)
    return extend_basis_to_target(S0, target_basis, l * m)


# ---------------------------------------------------------------------------
# Logical operators
# ---------------------------------------------------------------------------


def kernel_supported_ops_ifft(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]):
    """Return the old kernel-supported representative Z/X operators.

    These operators commute with the checks, but they are not guaranteed to be
    a complete or canonical logical basis modulo stabilizers.
    """
    S = kernel_basis_ifft(l, m, A_exp, B_exp, do_half_period_shifts=True)
    _, _, HX, HZ, C = build_blocks(l, m, A_exp, B_exp)
    Z_ops = np.hstack([S, np.zeros_like(S, dtype=np.uint8)])
    X_ops = np.hstack([np.zeros_like(S, dtype=np.uint8), gf2_matmul(S, C)])
    return Z_ops.astype(np.uint8), X_ops.astype(np.uint8), HX, HZ


def logical_ops_from_checks(HX: np.ndarray, HZ: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Construct a canonical CSS logical basis.

    Returns ``(Z_ops, X_ops)`` with one row per encoded qubit and with
    ``X_ops @ Z_ops.T == I`` over GF(2). Stabilizer rows are quotiented out.
    """
    HX = (np.asarray(HX, dtype=np.uint8) & 1)
    HZ = (np.asarray(HZ, dtype=np.uint8) & 1)
    if HX.ndim != 2 or HZ.ndim != 2 or HX.shape[1] != HZ.shape[1]:
        raise ValueError("HX and HZ must be binary matrices with the same number of columns.")
    n = HX.shape[1]

    # Optional consistency check. It is cheap for the code sizes used here.
    if gf2_matmul(HX, HZ.T).any():
        raise ValueError("HX and HZ do not commute: HX @ HZ.T is nonzero over GF(2).")

    k = int(n - gf2_rank(HX) - gf2_rank(HZ))
    if k < 0:
        raise ValueError("Invalid CSS check matrices: negative encoded dimension.")

    Z_kernel = gf2_nullspace(HX)
    X_kernel = gf2_nullspace(HZ)
    Z_stabs = gf2_rref_basis(HZ, n)
    X_stabs = gf2_rref_basis(HX, n)

    Z_ops = _quotient_basis(Z_kernel, Z_stabs, n)
    X_ops_raw = _quotient_basis(X_kernel, X_stabs, n)

    if Z_ops.shape[0] != k or X_ops_raw.shape[0] != k:
        raise RuntimeError(
            f"Failed to extract logical quotient bases: expected {k}, "
            f"got Z={Z_ops.shape[0]}, X={X_ops_raw.shape[0]}."
        )
    if k == 0:
        return Z_ops, X_ops_raw

    pairing = gf2_matmul(X_ops_raw, Z_ops.T)
    inv_pairing = gf2_inverse_square(pairing)
    X_ops = gf2_matmul(inv_pairing, X_ops_raw)

    final_pairing = gf2_matmul(X_ops, Z_ops.T)
    if not np.array_equal(final_pairing, np.eye(k, dtype=np.uint8)):
        raise RuntimeError("Internal error: logical pairing was not canonicalized.")
    return Z_ops.astype(np.uint8), X_ops.astype(np.uint8)


def logical_ops_canonical(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]):
    """Construct a full canonical logical Z/X basis for the BB code."""
    _, _, HX, HZ, _ = build_blocks(l, m, A_exp, B_exp)
    Z_ops, X_ops = logical_ops_from_checks(HX, HZ)
    return Z_ops, X_ops, HX, HZ


def logical_ops_ifft(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]):
    """Backward-compatible public name for a full canonical logical basis."""
    return logical_ops_canonical(l, m, A_exp, B_exp)


def validate_logical_basis(Z_ops: np.ndarray, X_ops: np.ndarray, HX: np.ndarray, HZ: np.ndarray) -> Dict[str, object]:
    """Validate CSS logical-basis conditions and return diagnostic data."""
    Z_ops = (np.asarray(Z_ops, dtype=np.uint8) & 1)
    X_ops = (np.asarray(X_ops, dtype=np.uint8) & 1)
    HX = (np.asarray(HX, dtype=np.uint8) & 1)
    HZ = (np.asarray(HZ, dtype=np.uint8) & 1)
    pairing = gf2_matmul(X_ops, Z_ops.T)
    eye = np.eye(Z_ops.shape[0], dtype=np.uint8) if Z_ops.shape[0] == X_ops.shape[0] else None
    return {
        "z_commutes_with_x_checks": not gf2_matmul(HX, Z_ops.T).any(),
        "x_commutes_with_z_checks": not gf2_matmul(HZ, X_ops.T).any(),
        "canonical_pairing": eye is not None and np.array_equal(pairing, eye),
        "pairing": pairing,
        "num_logicals": int(Z_ops.shape[0]),
        "blocklength": int(HX.shape[1]),
    }


# ---------------------------------------------------------------------------
# Example / smoke test
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    l_val, m_val = 21, 18
    A_exponents = [(3, 0), (0, 10), (0, 17)]
    B_exponents = [(0, 5), (3, 0), (19, 0)]

    k_fourier = k_bb_code_final(l_val, m_val, A_exponents, B_exponents)
    k_direct = direct_k_bb_code(l_val, m_val, A_exponents, B_exponents)
    print("k_fourier =", k_fourier)
    print("k_direct  =", k_direct)

    Zs, Xs, HX, HZ = logical_ops_canonical(l_val, m_val, A_exponents, B_exponents)
    report = validate_logical_basis(Zs, Xs, HX, HZ)
    print("Z_ops", Zs.shape, "X_ops", Xs.shape)
    print("valid =", {k: v for k, v in report.items() if k != "pairing"})
    assert k_fourier == k_direct == Zs.shape[0] == Xs.shape[0]
    assert report["z_commutes_with_x_checks"]
    assert report["x_commutes_with_z_checks"]
    assert report["canonical_pairing"]
