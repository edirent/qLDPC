"""
Logical-operator utilities for bivariate-bicycle (BB) CSS codes.

This module complements ``bb_fourier_tools_fast.py``.  It is not optimized for
speed; it focuses on producing a complete, paired, checkable logical-qubit
basis.

Main outputs
------------
* CSS supports: ``X_ops`` and ``Z_ops`` with shape (k, n).
* qLDPC-style symplectic matrix with shape (2*k, 2*n): first k rows are
  logical X operators, second k rows are logical Z operators.
* Per-logical-qubit support lists in the BB L/R torus coordinate convention.

Conventions
-----------
The physical qubit index range is 0..2*l*m-1.  The first l*m entries are L
qubits and the second l*m entries are R qubits.  A local index p*m+q means the
unit cell x^p y^q on the l-by-m torus.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Tuple
import json
import random

import numpy as np

from bb_fourier_tools_fast import (
    build_blocks,
    kernel_supported_ops_ifft,
    k_bb_code_final,
    normalize_exponents,
    gf2_matmul,
    gf2_rank,
    gf2_nullspace,
    gf2_rref_basis,
    gf2_inverse_square,
    _GF2RowSpace,
    _as_binary_matrix,
)

Exponent = Tuple[int, int]
PauliKind = Literal["X", "Z"]
PairingStrategy = Literal["auto", "transform_x", "transform_z"]


@dataclass(frozen=True)
class QubitSupport:
    """A physical-qubit support entry in BB torus coordinates."""

    index: int
    side: Literal["L", "R"]
    x: int
    y: int

    @property
    def monomial(self) -> str:
        """Return the local torus coordinate as a compact monomial label."""
        return _monomial_to_string((self.x, self.y))

    @property
    def label(self) -> str:
        """Return the side-qualified support label, e.g. ``L:x^2y``."""
        return f"{self.side}:{self.monomial}"


@dataclass(frozen=True)
class LogicalQubit:
    """One paired logical qubit, represented by one X and one Z support."""

    index: int
    x_support: Tuple[QubitSupport, ...]
    z_support: Tuple[QubitSupport, ...]

    @property
    def x_weight(self) -> int:
        """Hamming weight of the logical X representative."""
        return len(self.x_support)

    @property
    def z_weight(self) -> int:
        """Hamming weight of the logical Z representative."""
        return len(self.z_support)

    def to_dict(self) -> Dict[str, object]:
        """Serialize the logical pair and its support entries for JSON export."""
        return {
            "logical_index": self.index,
            "x_weight": self.x_weight,
            "z_weight": self.z_weight,
            "x_support": [entry.__dict__ | {"label": entry.label} for entry in self.x_support],
            "z_support": [entry.__dict__ | {"label": entry.label} for entry in self.z_support],
        }


def _monomial_to_string(exp: Exponent) -> str:
    """Format one normalized exponent pair as a human-readable monomial."""
    i, j = exp
    pieces: List[str] = []
    if i:
        pieces.append("x" if i == 1 else f"x^{i}")
    if j:
        pieces.append("y" if j == 1 else f"y^{j}")
    return "".join(pieces) if pieces else "1"


def _as_bin(M: Iterable, ncols: Optional[int] = None) -> np.ndarray:
    """Normalize matrix-like input to a binary uint8 NumPy matrix."""
    return _as_binary_matrix(M, ncols)


def _filter_kernel_rows(rows: np.ndarray, check: np.ndarray, ncols: int) -> np.ndarray:
    """Keep rows that commute with the relevant CSS checks."""
    rows = _as_bin(rows, ncols)
    if rows.shape[0] == 0:
        return rows
    check = _as_bin(check, ncols)
    ok = gf2_matmul(check, rows.T).sum(axis=0) == 0
    return rows[ok]


def _quotient_basis_seeded(
    space_basis: np.ndarray,
    mod_basis: np.ndarray,
    seeds: Optional[np.ndarray],
    ncols: int,
    target_dim: int,
) -> np.ndarray:
    """Basis of span(space_basis)/span(mod_basis), preferring seed rows.

    Rows are returned as representatives in the ambient ncols-dimensional space.
    The ordering is intentionally not RREF-normalized: keeping early seed rows
    gives more interpretable representatives.
    """
    rowspace = _GF2RowSpace(ncols)
    for row in _as_bin(mod_basis, ncols):
        rowspace.add(row)

    chosen: List[np.ndarray] = []
    candidates: List[np.ndarray] = []
    if seeds is not None:
        candidates.extend(list(_as_bin(seeds, ncols)))
    candidates.extend(list(_as_bin(space_basis, ncols)))

    for row in candidates:
        if len(chosen) >= target_dim:
            break
        before = rowspace.rank
        if rowspace.add(row) and rowspace.rank > before:
            chosen.append((row.copy() & 1).astype(np.uint8))

    if len(chosen) != target_dim:
        raise RuntimeError(f"Could not build a quotient basis of dimension {target_dim}; got {len(chosen)}.")
    return np.vstack(chosen).astype(np.uint8) if chosen else np.zeros((0, ncols), dtype=np.uint8)


def _candidate_generators(stabilizers: np.ndarray, ncols: int) -> np.ndarray:
    """Low-weight plus row-reduced stabilizer generators for coset reduction."""
    stabs = _as_bin(stabilizers, ncols)
    if stabs.shape[0] == 0:
        return stabs
    rref = gf2_rref_basis(stabs, ncols)
    gens = np.vstack([stabs, rref]).astype(np.uint8)
    # Remove duplicate rows while preserving a useful low-weight order.
    weights = gens.sum(axis=1)
    order = np.argsort(weights, kind="stable")
    seen = set()
    unique: List[np.ndarray] = []
    for idx in order:
        row = gens[int(idx)]
        key = row.tobytes()
        if key not in seen and row.any():
            seen.add(key)
            unique.append(row.copy())
    return np.vstack(unique).astype(np.uint8) if unique else np.zeros((0, ncols), dtype=np.uint8)


def reduce_by_stabilizers(
    op: np.ndarray,
    stabilizers: np.ndarray,
    *,
    trials: int = 64,
    seed: int = 0,
) -> np.ndarray:
    """Heuristically reduce the weight of one logical representative.

    This only adds stabilizer rows, so the logical class and pairings are
    preserved.  The problem of finding the true minimum-weight coset
    representative is generally hard; this function is a deterministic+random
    greedy heuristic intended to improve readability, not to certify distance.
    """
    v0 = (np.asarray(op, dtype=np.uint8).copy() & 1)
    ncols = int(v0.size)
    gens = _candidate_generators(stabilizers, ncols)
    if gens.shape[0] == 0:
        return v0

    def greedy(start: np.ndarray, order: Sequence[int]) -> np.ndarray:
        """Run one deterministic greedy descent over a chosen generator order."""
        v = start.copy()
        improved = True
        while improved:
            improved = False
            wt = int(v.sum())
            for idx in order:
                cand = v ^ gens[int(idx)]
                cand_wt = int(cand.sum())
                if cand_wt < wt:
                    v = cand
                    wt = cand_wt
                    improved = True
        return v

    base_order = list(range(gens.shape[0]))
    best = greedy(v0, base_order)
    best_wt = int(best.sum())

    rng = random.Random(seed)
    for _ in range(max(0, int(trials))):
        # Random start in the same coset, followed by greedy descent.
        start = v0.copy()
        # Sparse random perturbation avoids making every trial too noisy.
        for idx in base_order:
            if rng.random() < min(0.06, 8.0 / max(1, gens.shape[0])):
                start ^= gens[idx]
        order = base_order.copy()
        rng.shuffle(order)
        cand = greedy(start, order)
        cand_wt = int(cand.sum())
        if cand_wt < best_wt:
            best = cand
            best_wt = cand_wt
    return best.astype(np.uint8)


def reduce_ops_by_stabilizers(
    ops: np.ndarray,
    stabilizers: np.ndarray,
    *,
    trials: int = 64,
    seed: int = 0,
) -> np.ndarray:
    """Reduce every logical row independently by adding stabilizer generators."""
    ops = _as_bin(ops, None)
    out = []
    for i, row in enumerate(ops):
        out.append(reduce_by_stabilizers(row, stabilizers, trials=trials, seed=seed + 1009 * i))
    return np.vstack(out).astype(np.uint8) if out else ops


def logical_ops_from_checks_paired(
    HX: np.ndarray,
    HZ: np.ndarray,
    *,
    z_seeds: Optional[np.ndarray] = None,
    x_seeds: Optional[np.ndarray] = None,
    pairing_strategy: PairingStrategy = "auto",
    reduce: bool = True,
    reduction_trials: int = 64,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Construct paired CSS logicals from check matrices.

    Returns ``(Z_ops, X_ops)`` with one row per logical qubit and
    ``X_ops @ Z_ops.T == I`` over GF(2).  ``z_seeds`` and ``x_seeds`` can be
    used to prefer existing representatives, such as the old IFFT kernel
    representatives.
    """
    HX = _as_bin(HX, None)
    HZ = _as_bin(HZ, None)
    if HX.ndim != 2 or HZ.ndim != 2 or HX.shape[1] != HZ.shape[1]:
        raise ValueError("HX and HZ must be binary matrices with the same number of columns.")
    n = int(HX.shape[1])
    if gf2_matmul(HX, HZ.T).any():
        raise ValueError("HX and HZ do not commute over GF(2).")

    rank_x = gf2_rank(HX)
    rank_z = gf2_rank(HZ)
    k = int(n - rank_x - rank_z)
    if k < 0:
        raise ValueError("Invalid CSS code: negative encoded dimension.")
    if k == 0:
        empty = np.zeros((0, n), dtype=np.uint8)
        return empty, empty

    z_stabs = gf2_rref_basis(HZ, n)
    x_stabs = gf2_rref_basis(HX, n)
    z_kernel = gf2_nullspace(HX)
    x_kernel = gf2_nullspace(HZ)

    z_seed_rows = _filter_kernel_rows(z_seeds, HX, n) if z_seeds is not None else None
    x_seed_rows = _filter_kernel_rows(x_seeds, HZ, n) if x_seeds is not None else None

    Z_raw = _quotient_basis_seeded(z_kernel, z_stabs, z_seed_rows, n, k)
    X_raw = _quotient_basis_seeded(x_kernel, x_stabs, x_seed_rows, n, k)

    P = gf2_matmul(X_raw, Z_raw.T)
    invP = gf2_inverse_square(P)

    candidates: List[Tuple[str, np.ndarray, np.ndarray]] = []
    if pairing_strategy in ("auto", "transform_x"):
        candidates.append(("transform_x", Z_raw.copy(), gf2_matmul(invP, X_raw)))
    if pairing_strategy in ("auto", "transform_z"):
        candidates.append(("transform_z", gf2_matmul(invP.T, Z_raw), X_raw.copy()))

    best: Optional[Tuple[int, int, str, np.ndarray, np.ndarray]] = None
    for name, Z, X in candidates:
        if reduce:
            Z = reduce_ops_by_stabilizers(Z, HZ, trials=reduction_trials, seed=seed + 17)
            X = reduce_ops_by_stabilizers(X, HX, trials=reduction_trials, seed=seed + 31)
        report = validate_css_logicals(HX, HZ, X, Z)
        if not report["valid"]:
            raise RuntimeError(f"Candidate {name} failed validation: {report}")
        total_wt = int(Z.sum() + X.sum())
        max_wt = int(max(Z.sum(axis=1).max(initial=0), X.sum(axis=1).max(initial=0)))
        item = (total_wt, max_wt, name, Z, X)
        if best is None or item[:2] < best[:2]:
            best = item

    assert best is not None
    _, _, _, Z_best, X_best = best
    return Z_best.astype(np.uint8), X_best.astype(np.uint8)


def logical_ops_bb_paired(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    *,
    use_ifft_seeds: bool = False,
    pairing_strategy: PairingStrategy = "auto",
    reduce: bool = True,
    reduction_trials: int = 64,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct paired CSS logicals for a BB code.

    Returns ``(Z_ops, X_ops, HX, HZ)``.  The number of rows in ``Z_ops`` and
    ``X_ops`` is exactly the quantum dimension k.
    """
    _, _, HX, HZ, _ = build_blocks(l, m, A_exp, B_exp)
    z_seeds = x_seeds = None
    if use_ifft_seeds:
        try:
            z_seeds, x_seeds, _, _ = kernel_supported_ops_ifft(l, m, A_exp, B_exp)
        except Exception:
            # Seeds are optional.  Fallback remains exact, only less aligned with the old method.
            z_seeds = x_seeds = None
    Z_ops, X_ops = logical_ops_from_checks_paired(
        HX,
        HZ,
        z_seeds=z_seeds,
        x_seeds=x_seeds,
        pairing_strategy=pairing_strategy,
        reduce=reduce,
        reduction_trials=reduction_trials,
        seed=seed,
    )
    return Z_ops, X_ops, HX, HZ


def validate_css_logicals(
    HX: np.ndarray,
    HZ: np.ndarray,
    X_ops: np.ndarray,
    Z_ops: np.ndarray,
) -> Dict[str, object]:
    """Validate a paired CSS logical basis."""
    HX = _as_bin(HX, None)
    HZ = _as_bin(HZ, None)
    X_ops = _as_bin(X_ops, HX.shape[1])
    Z_ops = _as_bin(Z_ops, HX.shape[1])
    k = int(X_ops.shape[0])
    rank_hx = gf2_rank(HX)
    rank_hz = gf2_rank(HZ)
    expected_k = int(HX.shape[1] - rank_hx - rank_hz)
    pairing = gf2_matmul(X_ops, Z_ops.T)
    eye = np.eye(k, dtype=np.uint8) if X_ops.shape[0] == Z_ops.shape[0] else None
    z_synd = gf2_matmul(HX, Z_ops.T)
    x_synd = gf2_matmul(HZ, X_ops.T)
    z_new_rank = gf2_rank(np.vstack([HZ, Z_ops])) - rank_hz
    x_new_rank = gf2_rank(np.vstack([HX, X_ops])) - rank_hx
    result = {
        "valid": False,
        "expected_k": expected_k,
        "num_x": int(X_ops.shape[0]),
        "num_z": int(Z_ops.shape[0]),
        "z_commutes_with_x_checks": not z_synd.any(),
        "x_commutes_with_z_checks": not x_synd.any(),
        "canonical_pairing": eye is not None and np.array_equal(pairing, eye),
        "z_nontrivial_count": int(z_new_rank),
        "x_nontrivial_count": int(x_new_rank),
        "pairing": pairing,
        "x_weights": [int(w) for w in X_ops.sum(axis=1)],
        "z_weights": [int(w) for w in Z_ops.sum(axis=1)],
    }
    result["valid"] = (
        result["num_x"] == result["num_z"] == expected_k
        and result["z_commutes_with_x_checks"]
        and result["x_commutes_with_z_checks"]
        and result["canonical_pairing"]
        and result["z_nontrivial_count"] == expected_k
        and result["x_nontrivial_count"] == expected_k
    )
    return result


def qubit_support(index: int, l: int, m: int) -> QubitSupport:
    """Convert a flat physical-qubit index to BB L/R torus coordinates."""
    n0 = int(l) * int(m)
    index = int(index)
    if index < 0 or index >= 2 * n0:
        raise ValueError(f"Qubit index {index} outside 0..{2*n0-1}.")
    side: Literal["L", "R"] = "L" if index < n0 else "R"
    local = index if side == "L" else index - n0
    return QubitSupport(index=index, side=side, x=local // int(m), y=local % int(m))


def support_entries(op: np.ndarray, l: int, m: int) -> Tuple[QubitSupport, ...]:
    """Return coordinate support entries for the nonzero positions of ``op``."""
    op = (np.asarray(op, dtype=np.uint8).reshape(-1) & 1)
    return tuple(qubit_support(int(idx), l, m) for idx in np.flatnonzero(op))


def logical_qubits_from_ops(X_ops: np.ndarray, Z_ops: np.ndarray, l: int, m: int) -> Tuple[LogicalQubit, ...]:
    """Pair X/Z operator rows and attach BB coordinate support metadata."""
    X_ops = _as_bin(X_ops, None)
    Z_ops = _as_bin(Z_ops, X_ops.shape[1])
    if X_ops.shape != Z_ops.shape:
        raise ValueError("X_ops and Z_ops must have the same shape.")
    return tuple(
        LogicalQubit(
            index=i + 1,
            x_support=support_entries(X_ops[i], l, m),
            z_support=support_entries(Z_ops[i], l, m),
        )
        for i in range(X_ops.shape[0])
    )


def qldpc_style_logical_matrix(X_ops: np.ndarray, Z_ops: np.ndarray) -> np.ndarray:
    """Return a qLDPC/CSS-style symplectic logical matrix of shape (2k, 2n)."""
    X_ops = _as_bin(X_ops, None)
    Z_ops = _as_bin(Z_ops, X_ops.shape[1])
    if X_ops.shape != Z_ops.shape:
        raise ValueError("X_ops and Z_ops must have the same shape.")
    k, n = X_ops.shape
    zero = np.zeros((k, n), dtype=np.uint8)
    return np.vstack([np.hstack([X_ops, zero]), np.hstack([zero, Z_ops])]).astype(np.uint8)


def export_logical_qubits_json(logicals: Sequence[LogicalQubit], path: str) -> None:
    """Write paired logical-qubit support data as pretty-printed JSON."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump([q.to_dict() for q in logicals], f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Paper-specific gross/two-gross logical bases from Tour de gross Appendix A.1
# ---------------------------------------------------------------------------


def _poly_normalize(poly: Iterable[Sequence[int]], l: int, m: int) -> List[Exponent]:
    """Normalize BB polynomial exponents with GF(2) duplicate cancellation."""
    return normalize_exponents(poly, l, m)


def _poly_shift(poly: Iterable[Sequence[int]], shift: Sequence[int], l: int, m: int) -> List[Exponent]:
    """Translate every monomial in a BB polynomial by one torus shift."""
    sx, sy = int(shift[0]), int(shift[1])
    return _poly_normalize([(i + sx, j + sy) for i, j in poly], l, m)


def _poly_transpose(poly: Iterable[Sequence[int]], l: int, m: int) -> List[Exponent]:
    """Return the transpose polynomial obtained by inverting all exponents."""
    return _poly_normalize([(-int(i), -int(j)) for i, j in poly], l, m)


def _op_from_lr(l: int, m: int, left_poly: Iterable[Sequence[int]], right_poly: Iterable[Sequence[int]]) -> np.ndarray:
    """Build a length-``2*l*m`` support vector from left/right polynomial supports."""
    n0 = int(l) * int(m)
    out = np.zeros(2 * n0, dtype=np.uint8)
    for i, j in _poly_normalize(left_poly, l, m):
        out[i * int(m) + j] ^= 1
    for i, j in _poly_normalize(right_poly, l, m):
        out[n0 + i * int(m) + j] ^= 1
    return out


def _paper_polys(code: Literal["gross", "two-gross"]):
    """Return Appendix-A.1 polynomial data for a named Tour-de-gross code."""
    if code == "gross":
        p = [(4, 0), (5, 0), (6, 1), (4, 2), (5, 4), (6, 5)]
        q = [(3, 0), (4, 0), (3, 1), (3, 2), (4, 2), (3, 5)]
        r = [(0, 0), (8, 0), (1, 1), (9, 1), (3, 4), (11, 4)]
        s = [(1, 0), (9, 0), (4, 4), (8, 4), (0, 5), (8, 5)]
        alpha = [(0, 0), (3, 5), (11, 5), (10, 1), (5, 4), (4, 2)]
        beta = [(0, 0), (2, 4), (1, 2), (2, 5), (1, 1), (3, 1)]
    elif code == "two-gross":
        p = [(2, 0), (2, 2), (8, 2), (9, 2), (3, 3), (4, 3), (7, 4), (8, 6), (6, 7), (7, 11)]
        q = [(3, 2), (5, 3), (7, 3), (11, 3), (8, 4), (8, 5), (6, 7), (4, 8), (1, 9), (1, 10)]
        r = [(4, 2), (11, 2), (0, 5), (1, 5), (5, 5), (6, 5), (1, 8), (8, 8), (2, 11), (10, 11)]
        s = [(2, 6), (11, 6), (2, 9), (11, 9), (2, 10), (11, 10), (8, 7), (11, 7), (5, 8), (11, 8)]
        alpha = [(0, 0), (0, 3), (3, 7), (11, 11), (2, 9), (7, 4)]
        beta = [(0, 0), (1, 1), (4, 0), (5, 4), (4, 3), (3, 5)]
    else:
        raise ValueError("code must be 'gross' or 'two-gross'.")
    return p, q, r, s, alpha, beta


def tour_de_gross_logicals(
    code: Literal["gross", "two-gross"],
    *,
    reduce: bool = False,
    reduction_trials: int = 64,
    seed: int = 0,
) -> Tuple[int, int, List[Exponent], List[Exponent], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the Appendix-A.1 gross/two-gross logical basis.

    This basis uses the Tour-de-gross polynomial convention
    A = 1 + y + x^3 y^-1, B = 1 + x + x^-1 y^-3.

    Returns ``(l, m, A_exp, B_exp, Z_ops, X_ops, HX, HZ)``.  The row order is
    logical qubits 1..12 from Eq. (34); ``X_ops[i]`` pairs with ``Z_ops[i]``.
    """
    if code == "gross":
        l, m = 12, 6
    elif code == "two-gross":
        l, m = 12, 12
    else:
        raise ValueError("code must be 'gross' or 'two-gross'.")

    A_exp = _poly_normalize([(0, 0), (0, 1), (3, -1)], l, m)
    B_exp = _poly_normalize([(0, 0), (1, 0), (-1, -3)], l, m)
    p, q, r, s, alpha, beta = _paper_polys(code)
    mu = nu = (1, 1)

    X_rows: List[np.ndarray] = []
    Z_rows: List[np.ndarray] = []

    qT = _poly_transpose(q, l, m)
    pT = _poly_transpose(p, l, m)
    sT = _poly_transpose(s, l, m)
    rT = _poly_transpose(r, l, m)

    for a, b in zip(alpha, beta):
        # X_i = X(alpha_i p, alpha_i q)
        X_rows.append(_op_from_lr(l, m, _poly_shift(p, a, l, m), _poly_shift(q, a, l, m)))
        # Z_i = Z(beta_i nu s^T, beta_i nu r^T)
        bnu = (b[0] + nu[0], b[1] + nu[1])
        Z_rows.append(_op_from_lr(l, m, _poly_shift(sT, bnu, l, m), _poly_shift(rT, bnu, l, m)))

    for a, b in zip(alpha, beta):
        bT = (-b[0], -b[1])
        aTmu = (-a[0] + mu[0], -a[1] + mu[1])
        # X_{i+6} = X(beta_i^T r, beta_i^T s)
        X_rows.append(_op_from_lr(l, m, _poly_shift(r, bT, l, m), _poly_shift(s, bT, l, m)))
        # Z_{i+6} = Z(alpha_i^T mu q^T, alpha_i^T mu p^T)
        Z_rows.append(_op_from_lr(l, m, _poly_shift(qT, aTmu, l, m), _poly_shift(pT, aTmu, l, m)))

    X_ops = np.vstack(X_rows).astype(np.uint8)
    Z_ops = np.vstack(Z_rows).astype(np.uint8)
    _, _, HX, HZ, _ = build_blocks(l, m, A_exp, B_exp)

    if reduce:
        # Normally disabled: the paper basis is intentionally structured for LPU use.
        X_ops = reduce_ops_by_stabilizers(X_ops, HX, trials=reduction_trials, seed=seed + 1)
        Z_ops = reduce_ops_by_stabilizers(Z_ops, HZ, trials=reduction_trials, seed=seed + 2)

    report = validate_css_logicals(HX, HZ, X_ops, Z_ops)
    if not report["valid"]:
        raise RuntimeError(f"Tour-de-gross {code} basis failed validation: {report}")
    return l, m, A_exp, B_exp, Z_ops, X_ops, HX, HZ


def summarize_logicals(X_ops: np.ndarray, Z_ops: np.ndarray) -> Dict[str, object]:
    """Summarize dimensions and support-weight ranges for a logical basis."""
    X_ops = _as_bin(X_ops, None)
    Z_ops = _as_bin(Z_ops, X_ops.shape[1])
    xw = X_ops.sum(axis=1).astype(int) if X_ops.shape[0] else np.array([], dtype=int)
    zw = Z_ops.sum(axis=1).astype(int) if Z_ops.shape[0] else np.array([], dtype=int)
    return {
        "k": int(X_ops.shape[0]),
        "n": int(X_ops.shape[1]),
        "x_weights": [int(w) for w in xw],
        "z_weights": [int(w) for w in zw],
        "x_weight_min": int(xw.min()) if xw.size else 0,
        "x_weight_max": int(xw.max()) if xw.size else 0,
        "z_weight_min": int(zw.min()) if zw.size else 0,
        "z_weight_max": int(zw.max()) if zw.size else 0,
        "total_weight": int(X_ops.sum() + Z_ops.sum()),
    }

# ---------------------------------------------------------------------------
# BB-specific upgrade of the old kernel-supported construction
# ---------------------------------------------------------------------------


def kernel_symmetric_logical_basis(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    *,
    reduce_weight: bool = True,
) -> object:
    """Compatibility wrapper for the maintained cleaned constructor.

    The experimental kernel-symmetric code path originally kept here depended
    on helper names that moved into ``bb_logical_tools_clean.py``.  The cleaned
    canonical constructor is now the validated implementation, so this legacy
    entry point delegates there while preserving the public call signature.
    """
    import bb_logical_tools_clean as clean

    return clean.canonical_logical_basis(
        l,
        m,
        A_exp,
        B_exp,
        use_ifft_seeds=True,
        reduce_weight=reduce_weight,
    )


def logical_basis_best_effort(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    *,
    prefer_tour_named_basis: bool = True,
    prefer_kernel_symmetric: bool = True,
    reduce_weight: bool = True,
) -> object:
    """Best-effort logical basis selector kept for backward compatibility.

    The maintained selector lives in ``bb_logical_tools_clean.py``.  This
    wrapper keeps older imports working, maps the Tour-de-gross preference, and
    leaves ``prefer_kernel_symmetric`` as a no-op compatibility flag.
    """
    import bb_logical_tools_clean as clean

    if prefer_tour_named_basis:
        return clean.best_logical_basis(
            l,
            m,
            A_exp,
            B_exp,
            prefer_tour_basis=True,
            reduce_weight=reduce_weight,
        )
    return clean.canonical_logical_basis(
        l,
        m,
        A_exp,
        B_exp,
        reduce_weight=reduce_weight,
    )
