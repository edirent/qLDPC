"""Clean logical-operator utilities for bivariate-bicycle (BB) CSS codes.

This file is correctness/readability oriented.  It complements the fast Fourier
``k`` routines by constructing *paired logical qubits*:

    X_ops[i] @ Z_ops[j].T = delta_ij  over GF(2)

where each row has length n = 2*l*m and stores the physical support on the BB
L/R qubit blocks.  A qLDPC-style symplectic matrix can be exported as
[X_support | Z_support] rows: first all logical X_i, then all logical Z_i.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Tuple
import argparse
import json
import os
import random
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import bb_fourier_tools_fast as bb

Exponent = Tuple[int, int]


# ---------------------------------------------------------------------------
# GF(2) helpers
# ---------------------------------------------------------------------------


def _bin(M: np.ndarray | Sequence[Sequence[int]] | Sequence[int]) -> np.ndarray:
    """Return a uint8 NumPy array with every entry reduced modulo two."""
    return (np.asarray(M, dtype=np.uint8) & 1).astype(np.uint8)


def gf2_matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Matrix product over GF(2), delegated to the shared Fourier helper."""
    return bb.gf2_matmul(_bin(A), _bin(B))


def gf2_rank(M: np.ndarray) -> int:
    """Rank of a binary matrix over GF(2)."""
    return int(bb.gf2_rank(_bin(M)))


def gf2_nullspace(M: np.ndarray) -> np.ndarray:
    """Return a row basis for the binary nullspace of ``M``."""
    return bb.gf2_nullspace(_bin(M))


def gf2_rref_basis(M: np.ndarray, ncols: Optional[int] = None) -> np.ndarray:
    """Return an independent row basis in reduced row-echelon form."""
    M = _bin(M)
    if M.size == 0:
        if ncols is None:
            ncols = 0
        return np.zeros((0, int(ncols)), dtype=np.uint8)
    if M.ndim == 1:
        M = M.reshape(1, -1)
    if ncols is None:
        ncols = M.shape[1]
    return bb.gf2_rref_basis(M, int(ncols))


def gf2_inverse_square(M: np.ndarray) -> np.ndarray:
    """Invert a nonsingular square binary matrix."""
    return bb.gf2_inverse_square(_bin(M))


class GF2RowSpace:
    """Incrementally maintain row-space membership over GF(2).

    Rows are stored by their leading pivot.  Adding a row reduces it against the
    existing pivots; the method returns ``True`` only when the row increases the
    represented rank.
    """

    def __init__(self, ncols: int):
        """Create an empty row space for vectors of length ``ncols``."""
        self.ncols = int(ncols)
        self.rows_by_pivot: Dict[int, np.ndarray] = {}

    @property
    def rank(self) -> int:
        """Current number of independent rows."""
        return len(self.rows_by_pivot)

    def add(self, row: np.ndarray) -> bool:
        """Add ``row`` if it is independent of the stored pivot rows."""
        v = _bin(row).reshape(-1).copy()
        if v.size != self.ncols:
            raise ValueError(f"Expected row length {self.ncols}; got {v.size}")
        while True:
            nz = np.flatnonzero(v)
            if nz.size == 0:
                return False
            piv = int(nz[0])
            old = self.rows_by_pivot.get(piv)
            if old is None:
                self.rows_by_pivot[piv] = v
                return True
            v ^= old


def quotient_basis(space_basis: np.ndarray, mod_basis: np.ndarray, ncols: int, seeds: Optional[np.ndarray] = None) -> np.ndarray:
    """Return representatives for span(space_basis)/span(mod_basis)."""
    rowspace = GF2RowSpace(ncols)
    mod = gf2_rref_basis(mod_basis, ncols)
    for row in mod:
        rowspace.add(row)

    chosen: List[np.ndarray] = []
    candidates: List[np.ndarray] = []
    if seeds is not None and np.asarray(seeds).size:
        candidates.extend(list(_bin(seeds).reshape(-1, ncols)))
    if np.asarray(space_basis).size:
        candidates.extend(list(gf2_rref_basis(space_basis, ncols)))

    for row in candidates:
        before = rowspace.rank
        if rowspace.add(row) and rowspace.rank > before:
            chosen.append(row.copy())
    if not chosen:
        return np.zeros((0, ncols), dtype=np.uint8)
    return np.vstack(chosen).astype(np.uint8)


# ---------------------------------------------------------------------------
# Polynomial and support helpers
# ---------------------------------------------------------------------------


def normalize_terms(terms: Iterable[Sequence[int]], l: int, m: int) -> List[Exponent]:
    """Normalize exponents modulo x^l=1, y^m=1; duplicate terms cancel."""
    l, m = int(l), int(m)
    if l <= 0 or m <= 0:
        raise ValueError("l and m must be positive.")
    counts: Dict[Exponent, int] = {}
    for term in terms:
        if len(term) != 2:
            raise ValueError(f"Bad exponent pair: {term!r}")
        key = (int(term[0]) % l, int(term[1]) % m)
        counts[key] = counts.get(key, 0) ^ 1
    return sorted([key for key, val in counts.items() if val])


def transpose_terms(terms: Iterable[Sequence[int]], l: int, m: int) -> List[Exponent]:
    """Invert all monomial exponents, matching polynomial transpose over the torus."""
    return normalize_terms(((-int(i), -int(j)) for i, j in terms), l, m)


def shift_terms(terms: Iterable[Sequence[int]], shift: Sequence[int], l: int, m: int) -> List[Exponent]:
    """Translate a polynomial support by ``shift=(dx, dy)`` on the torus."""
    dx, dy = int(shift[0]), int(shift[1])
    return normalize_terms(((int(i) + dx, int(j) + dy) for i, j in terms), l, m)


def terms_to_vec(terms: Iterable[Sequence[int]], l: int, m: int) -> np.ndarray:
    """Encode normalized monomial support as a flattened ``l*m`` binary vector."""
    v = np.zeros(int(l) * int(m), dtype=np.uint8)
    for i, j in normalize_terms(terms, l, m):
        v[i * int(m) + j] ^= 1
    return v


def vec_to_terms(vec: np.ndarray, l: int, m: int) -> List[Exponent]:
    """Decode a flattened torus support vector back to exponent pairs."""
    vec = _bin(vec).reshape(int(l) * int(m))
    return [(int(idx // int(m)), int(idx % int(m))) for idx in np.flatnonzero(vec)]


def lr_vec(left_terms: Iterable[Sequence[int]], right_terms: Iterable[Sequence[int]], l: int, m: int) -> np.ndarray:
    """Concatenate left and right BB torus supports into one physical vector."""
    return np.hstack([terms_to_vec(left_terms, l, m), terms_to_vec(right_terms, l, m)]).astype(np.uint8)


def split_lr(vec: np.ndarray, l: int, m: int) -> Tuple[np.ndarray, np.ndarray]:
    """Split a physical support vector into its left and right torus halves."""
    n0 = int(l) * int(m)
    vec = _bin(vec).reshape(2 * n0)
    return vec[:n0], vec[n0:]


def _monomial_str(term: Exponent) -> str:
    """Format one exponent pair as ``1``, ``x^i``, ``y^j``, or ``x^iy^j``."""
    i, j = term
    pieces = []
    if i:
        pieces.append("x" if i == 1 else f"x^{i}")
    if j:
        pieces.append("y" if j == 1 else f"y^{j}")
    return "".join(pieces) if pieces else "1"


def terms_to_poly_str(terms: Iterable[Sequence[int]], l: int, m: int) -> str:
    """Format a polynomial support as a readable sum of monomials."""
    norm = normalize_terms(terms, l, m)
    return "0" if not norm else " + ".join(_monomial_str(t) for t in norm)


def lr_vec_to_poly_str(vec: np.ndarray, l: int, m: int) -> str:
    """Format a full left/right support vector as two polynomial strings."""
    left, right = split_lr(vec, l, m)
    return f"L: {terms_to_poly_str(vec_to_terms(left, l, m), l, m)} | R: {terms_to_poly_str(vec_to_terms(right, l, m), l, m)}"


def support_entries(vec: np.ndarray, l: int, m: int) -> List[Dict[str, object]]:
    """Return JSON-friendly support entries with side, index, and torus labels."""
    left, right = split_lr(vec, l, m)
    out: List[Dict[str, object]] = []
    for side, half, offset in [("L", left, 0), ("R", right, int(l) * int(m))]:
        for i, j in vec_to_terms(half, l, m):
            idx = offset + i * int(m) + j
            out.append({"index": idx, "side": side, "x": i, "y": j, "label": f"{side}:{_monomial_str((i,j))}"})
    return out


# ---------------------------------------------------------------------------
# Validation and container
# ---------------------------------------------------------------------------


@dataclass
class LogicalBasis:
    """Complete paired logical basis for one BB CSS code.

    ``X_ops`` and ``Z_ops`` have one row per logical qubit and length
    ``n = 2*l*m``.  The rows are paired so that ``X_ops @ Z_ops.T`` is the
    identity over GF(2), while ``HX`` and ``HZ`` retain the stabilizer checks
    used to validate commutation and quotient nontriviality.
    """

    l: int
    m: int
    A_exp: List[Exponent]
    B_exp: List[Exponent]
    Z_ops: np.ndarray  # shape (k,n), Z physical support
    X_ops: np.ndarray  # shape (k,n), X physical support
    HX: np.ndarray
    HZ: np.ndarray
    labels: List[str]
    source: str
    notes: str = ""

    @property
    def n(self) -> int:
        """Number of physical qubits in the BB code."""
        return int(self.HX.shape[1])

    @property
    def k(self) -> int:
        """Number of encoded logical qubits represented by this basis."""
        return int(self.X_ops.shape[0])

    @property
    def pairing(self) -> np.ndarray:
        """Binary symplectic X/Z pairing matrix for the stored representatives."""
        return gf2_matmul(self.X_ops, self.Z_ops.T)

    def validate(self) -> Dict[str, object]:
        """Run all CSS logical-basis consistency checks."""
        return validate_css_logicals(self.HX, self.HZ, self.X_ops, self.Z_ops)

    def qldpc_style_logical_matrix(self) -> np.ndarray:
        """Return shape (2k, 2n): [X|0] rows followed by [0|Z] rows."""
        k, n = self.X_ops.shape
        zero = np.zeros((k, n), dtype=np.uint8)
        return np.vstack([np.hstack([self.X_ops, zero]), np.hstack([zero, self.Z_ops])]).astype(np.uint8)

    def support_table(self, *, include_entries: bool = False) -> List[Dict[str, object]]:
        """Return one row of weights and polynomial supports per logical qubit.

        ``include_entries`` expands each polynomial support into individual
        physical-qubit entries, which is useful for JSON but verbose for reports.
        """
        rows: List[Dict[str, object]] = []
        for i, label in enumerate(self.labels):
            row = {
                "logical": label,
                "X_weight": int(self.X_ops[i].sum()),
                "Z_weight": int(self.Z_ops[i].sum()),
                "X_polynomial": lr_vec_to_poly_str(self.X_ops[i], self.l, self.m),
                "Z_polynomial": lr_vec_to_poly_str(self.Z_ops[i], self.l, self.m),
            }
            if include_entries:
                row["X_support"] = support_entries(self.X_ops[i], self.l, self.m)
                row["Z_support"] = support_entries(self.Z_ops[i], self.l, self.m)
            rows.append(row)
        return rows

    def to_json_dict(self, *, include_entries: bool = True) -> Dict[str, object]:
        """Serialize the basis, metadata, and validation report as plain objects."""
        report = self.validate().copy()
        if isinstance(report.get("pairing"), np.ndarray):
            report["pairing"] = report["pairing"].tolist()
        return {
            "l": self.l,
            "m": self.m,
            "n": self.n,
            "k": self.k,
            "A_exp": self.A_exp,
            "B_exp": self.B_exp,
            "source": self.source,
            "notes": self.notes,
            "validation": report,
            "logical_qubits": self.support_table(include_entries=include_entries),
        }

    def write_json(self, path: str, *, include_entries: bool = True) -> None:
        """Write a JSON report for the basis."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_json_dict(include_entries=include_entries), f, indent=2, ensure_ascii=False)

    def write_markdown(self, path: str, *, max_poly_chars: int = 240) -> None:
        """Write a compact Markdown table of logical weights and supports."""
        def trim(s: str) -> str:
            """Truncate long polynomial strings so Markdown tables stay readable."""
            return s if len(s) <= max_poly_chars else s[: max_poly_chars - 3] + "..."

        report = self.validate()
        lines = [
            f"# Logical basis report: {self.source}",
            "",
            f"- l={self.l}, m={self.m}, n={self.n}, k={self.k}",
            f"- A exponents: {self.A_exp}",
            f"- B exponents: {self.B_exp}",
            f"- notes: {self.notes}",
            "",
            "## Validation",
        ]
        for key, val in report.items():
            if key != "pairing":
                lines.append(f"- {key}: {val}")
        lines += ["", "## Logical qubits", "| logical | wt(X) | X support | wt(Z) | Z support |", "|---|---:|---|---:|---|"]
        for row in self.support_table(include_entries=False):
            lines.append(
                f"| {row['logical']} | {row['X_weight']} | `{trim(str(row['X_polynomial']))}` | "
                f"{row['Z_weight']} | `{trim(str(row['Z_polynomial']))}` |"
            )
        lines.append("")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))


def validate_css_logicals(HX: np.ndarray, HZ: np.ndarray, X_ops: np.ndarray, Z_ops: np.ndarray) -> Dict[str, object]:
    """Validate that X/Z rows form a complete canonical CSS logical basis.

    The report checks four independent properties: both Pauli types commute
    with the opposite stabilizer checks, the X/Z pairing matrix is identity,
    and each set contributes exactly ``k`` new quotient directions modulo its
    corresponding stabilizer space.
    """
    HX = _bin(HX)
    HZ = _bin(HZ)
    X_ops = _bin(X_ops)
    Z_ops = _bin(Z_ops)
    if X_ops.ndim != 2 or Z_ops.ndim != 2 or X_ops.shape != Z_ops.shape:
        raise ValueError("X_ops and Z_ops must have the same 2D shape.")
    if HX.shape[1] != X_ops.shape[1] or HZ.shape[1] != X_ops.shape[1]:
        raise ValueError("Check matrices and logical operators have incompatible widths.")
    rank_hx = gf2_rank(HX)
    rank_hz = gf2_rank(HZ)
    expected_k = int(HX.shape[1] - rank_hx - rank_hz)
    k = int(X_ops.shape[0])
    pairing = gf2_matmul(X_ops, Z_ops.T)
    z_synd = gf2_matmul(HX, Z_ops.T)
    x_synd = gf2_matmul(HZ, X_ops.T)
    z_inc = gf2_rank(np.vstack([HZ, Z_ops])) - rank_hz
    x_inc = gf2_rank(np.vstack([HX, X_ops])) - rank_hx
    result: Dict[str, object] = {
        "valid": False,
        "expected_k": expected_k,
        "num_x": int(X_ops.shape[0]),
        "num_z": int(Z_ops.shape[0]),
        "z_commutes_with_x_checks": not z_synd.any(),
        "x_commutes_with_z_checks": not x_synd.any(),
        "canonical_pairing": np.array_equal(pairing, np.eye(k, dtype=np.uint8)),
        "z_nonstabilizer_quotient_rank": int(z_inc),
        "x_nonstabilizer_quotient_rank": int(x_inc),
        "pairing_rank": gf2_rank(pairing),
        "pairing": pairing,
        "X_weights": [int(x.sum()) for x in X_ops],
        "Z_weights": [int(z.sum()) for z in Z_ops],
    }
    result["valid"] = bool(
        result["num_x"] == result["num_z"] == expected_k
        and result["z_commutes_with_x_checks"]
        and result["x_commutes_with_z_checks"]
        and result["canonical_pairing"]
        and result["z_nonstabilizer_quotient_rank"] == expected_k
        and result["x_nonstabilizer_quotient_rank"] == expected_k
    )
    return result


# ---------------------------------------------------------------------------
# Construction and optimization
# ---------------------------------------------------------------------------


def _kernel_seeds(l: int, m: int, A_exp: List[Exponent], B_exp: List[Exponent]) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Try to reuse IFFT-supported representatives as ordering seeds."""
    try:
        Z_old, X_old, _, _ = bb.kernel_supported_ops_ifft(l, m, A_exp, B_exp)
        return _bin(Z_old), _bin(X_old)
    except Exception:
        return None, None


def canonical_css_logicals_from_checks(
    HX: np.ndarray,
    HZ: np.ndarray,
    *,
    z_seeds: Optional[np.ndarray] = None,
    x_seeds: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (Z_ops, X_ops) with X_ops @ Z_ops.T = I."""
    HX = _bin(HX)
    HZ = _bin(HZ)
    if gf2_matmul(HX, HZ.T).any():
        raise ValueError("HX and HZ do not commute.")
    n = int(HX.shape[1])
    k = int(n - gf2_rank(HX) - gf2_rank(HZ))
    if k < 0:
        raise ValueError("Negative encoded dimension.")
    if k == 0:
        empty = np.zeros((0, n), dtype=np.uint8)
        return empty, empty

    Z_kernel = gf2_nullspace(HX)
    X_kernel = gf2_nullspace(HZ)
    Z_stabs = gf2_rref_basis(HZ, n)
    X_stabs = gf2_rref_basis(HX, n)
    Z_raw = quotient_basis(Z_kernel, Z_stabs, n, seeds=z_seeds)
    X_raw = quotient_basis(X_kernel, X_stabs, n, seeds=x_seeds)

    # Truncate to the quotient dimension in case seeds+candidates produced more.
    Z_raw = Z_raw[:k]
    X_raw = X_raw[:k]
    if Z_raw.shape[0] != k or X_raw.shape[0] != k:
        raise RuntimeError(f"Expected k={k}, got Z={Z_raw.shape[0]}, X={X_raw.shape[0]}.")

    P = gf2_matmul(X_raw, Z_raw.T)
    if gf2_rank(P) != k:
        raise RuntimeError("X/Z quotient pairing is singular.")

    # Compare transforming X versus transforming Z; choose lower total weight.
    invP = gf2_inverse_square(P)
    cand1 = (Z_raw.copy(), gf2_matmul(invP, X_raw))         # X <- P^-1 X
    cand2 = (gf2_matmul(invP.T, Z_raw), X_raw.copy())       # Z <- P^-T Z
    if int(cand2[0].sum() + cand2[1].sum()) < int(cand1[0].sum() + cand1[1].sum()):
        return cand2[0].astype(np.uint8), cand2[1].astype(np.uint8)
    return cand1[0].astype(np.uint8), cand1[1].astype(np.uint8)


def _candidate_stabilizers(stabs: np.ndarray, ncols: int) -> np.ndarray:
    """Prepare stabilizer rows for greedy coset-weight cleanup.

    Original stabilizer rows and their RREF basis are combined, sorted by
    weight, and deduplicated.  The resulting order makes the first greedy pass
    deterministic while still keeping useful low-weight generators early.
    """
    stabs = _bin(stabs).reshape(-1, ncols)
    if stabs.size == 0:
        return np.zeros((0, ncols), dtype=np.uint8)
    rref = gf2_rref_basis(stabs, ncols)
    gens = np.vstack([stabs, rref])
    order = np.argsort(gens.sum(axis=1), kind="stable")
    seen = set()
    rows: List[np.ndarray] = []
    for idx in order:
        row = gens[int(idx)]
        if not row.any():
            continue
        key = row.tobytes()
        if key not in seen:
            seen.add(key)
            rows.append(row.copy())
    return np.vstack(rows).astype(np.uint8) if rows else np.zeros((0, ncols), dtype=np.uint8)


def reduce_one_by_stabilizers(op: np.ndarray, stabs: np.ndarray, *, trials: int = 32, seed: int = 0) -> np.ndarray:
    """Heuristically lower one logical representative's support weight.

    Only stabilizer rows are added, so the logical coset is unchanged.  The
    randomized trials are reproducible and are a heuristic cleanup rather than a
    minimum-distance certificate.
    """
    v0 = _bin(op).reshape(-1)
    gens = _candidate_stabilizers(stabs, v0.size)
    if gens.shape[0] == 0:
        return v0.copy()

    def greedy(start: np.ndarray, order: Sequence[int]) -> np.ndarray:
        """Greedily add stabilizers whenever they strictly reduce weight."""
        v = start.copy()
        changed = True
        while changed:
            changed = False
            current = int(v.sum())
            for idx in order:
                cand = v ^ gens[int(idx)]
                wt = int(cand.sum())
                if wt < current:
                    v, current, changed = cand, wt, True
        return v

    base_order = list(range(gens.shape[0]))
    best = greedy(v0, base_order)
    best_wt = int(best.sum())
    rng = random.Random(seed)
    for _ in range(int(trials)):
        start = v0.copy()
        for idx in base_order:
            if rng.random() < min(0.06, 8.0 / max(1, gens.shape[0])):
                start ^= gens[idx]
        order = base_order.copy()
        rng.shuffle(order)
        cand = greedy(start, order)
        wt = int(cand.sum())
        if wt < best_wt:
            best, best_wt = cand, wt
    return best.astype(np.uint8)


def reduce_ops_by_stabilizers(ops: np.ndarray, stabs: np.ndarray, *, trials: int = 32, seed: int = 0) -> np.ndarray:
    """Apply stabilizer-coset reduction to every row of ``ops``."""
    ops = _bin(ops)
    if ops.shape[0] == 0:
        return ops
    return np.vstack([reduce_one_by_stabilizers(row, stabs, trials=trials, seed=seed + 1009*i) for i, row in enumerate(ops)]).astype(np.uint8)


def optimize_logicals_by_stabilizers(Z_ops: np.ndarray, X_ops: np.ndarray, HX: np.ndarray, HZ: np.ndarray, *, trials: int = 32) -> Tuple[np.ndarray, np.ndarray]:
    """Reduce X and Z representatives, then revalidate the paired basis."""
    X_red = reduce_ops_by_stabilizers(X_ops, HX, trials=trials, seed=11)
    Z_red = reduce_ops_by_stabilizers(Z_ops, HZ, trials=trials, seed=23)
    report = validate_css_logicals(HX, HZ, X_red, Z_red)
    if not report["valid"]:
        raise RuntimeError("Stabilizer reduction broke logical validation.")
    return Z_red, X_red


def canonical_logical_basis(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    *,
    use_ifft_seeds: bool = True,
    reduce_weight: bool = True,
) -> LogicalBasis:
    """Build a generic canonical logical basis from BB polynomial data.

    The construction forms the CSS check matrices, selects quotient-basis
    representatives from the two kernels, canonicalizes the X/Z pairing, and
    optionally reduces representatives by adding stabilizers.
    """
    l, m = int(l), int(m)
    A_norm = normalize_terms(A_exp, l, m)
    B_norm = normalize_terms(B_exp, l, m)
    _, _, HX, HZ, _ = bb.build_blocks(l, m, A_norm, B_norm)
    z_seed = x_seed = None
    if use_ifft_seeds:
        z_seed, x_seed = _kernel_seeds(l, m, A_norm, B_norm)
    Z_ops, X_ops = canonical_css_logicals_from_checks(HX, HZ, z_seeds=z_seed, x_seeds=x_seed)
    if reduce_weight:
        Z_ops, X_ops = optimize_logicals_by_stabilizers(Z_ops, X_ops, HX, HZ)
    basis = LogicalBasis(
        l=l, m=m, A_exp=A_norm, B_exp=B_norm, Z_ops=Z_ops, X_ops=X_ops, HX=HX, HZ=HZ,
        labels=[f"q{i}" for i in range(1, X_ops.shape[0] + 1)],
        source="canonical quotient basis" + (" with IFFT seeds" if use_ifft_seeds else "") + (" and stabilizer weight cleanup" if reduce_weight else ""),
        notes="Generic algebraic basis. Labels are canonical-basis labels, not architecture/LPU labels unless a paper basis is used.",
    )
    if not basis.validate()["valid"]:
        raise RuntimeError("Canonical logical basis failed validation.")
    return basis


# ---------------------------------------------------------------------------
# Old representatives diagnostic
# ---------------------------------------------------------------------------


def legacy_kernel_representatives(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the older IFFT-kernel representatives for comparison only."""
    return bb.kernel_supported_ops_ifft(l, m, A_exp, B_exp)


def diagnose_legacy_logicals(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> Dict[str, object]:
    """Measure why legacy representatives are not necessarily a full paired basis."""
    Z_old, X_old, HX, HZ = legacy_kernel_representatives(l, m, A_exp, B_exp)
    n = HX.shape[1]
    k = int(n - gf2_rank(HX) - gf2_rank(HZ))
    pairing = gf2_matmul(X_old, Z_old.T) if X_old.size and Z_old.size else np.zeros((X_old.shape[0], Z_old.shape[0]), dtype=np.uint8)
    return {
        "n": int(n),
        "expected_k": k,
        "old_Z_rows": int(Z_old.shape[0]),
        "old_X_rows": int(X_old.shape[0]),
        "old_rows_equal_k": bool(Z_old.shape[0] == X_old.shape[0] == k),
        "old_z_commutes_with_x_checks": not gf2_matmul(HX, Z_old.T).any(),
        "old_x_commutes_with_z_checks": not gf2_matmul(HZ, X_old.T).any(),
        "old_Z_rank_increment_mod_Z_stabilizers": int(gf2_rank(np.vstack([HZ, Z_old])) - gf2_rank(HZ)),
        "old_X_rank_increment_mod_X_stabilizers": int(gf2_rank(np.vstack([HX, X_old])) - gf2_rank(HX)),
        "old_pairing_rank": gf2_rank(pairing),
        "old_pairing": pairing,
    }


# ---------------------------------------------------------------------------
# Tour-de-gross explicit gross/two-gross bases
# ---------------------------------------------------------------------------


_TOUR_A = [(0, 0), (0, 1), (3, -1)]
_TOUR_B = [(0, 0), (1, 0), (-1, -3)]


def _tour_polys(code: Literal["gross", "two-gross"]):
    """Return the Appendix-A.1 data for the named Tour-de-gross construction."""
    if code == "gross":
        l, m = 12, 6
        p = [(4, 0), (5, 0), (6, 1), (4, 2), (5, 4), (6, 5)]
        q = [(3, 0), (4, 0), (3, 1), (3, 2), (4, 2), (3, 5)]
        r = [(0, 0), (8, 0), (1, 1), (9, 1), (3, 4), (11, 4)]
        s = [(1, 0), (9, 0), (4, 4), (8, 4), (0, 5), (8, 5)]
        alpha = [(0, 0), (3, 5), (11, 5), (10, 1), (5, 4), (4, 2)]
        beta = [(0, 0), (2, 4), (1, 2), (2, 5), (1, 1), (3, 1)]
    elif code == "two-gross":
        l, m = 12, 12
        p = [(2, 0), (2, 2), (8, 2), (9, 2), (3, 3), (4, 3), (7, 4), (8, 6), (6, 7), (7, 11)]
        q = [(3, 2), (5, 3), (7, 3), (11, 3), (8, 4), (8, 5), (6, 7), (4, 8), (1, 9), (1, 10)]
        r = [(4, 2), (11, 2), (0, 5), (1, 5), (5, 5), (6, 5), (1, 8), (8, 8), (2, 11), (10, 11)]
        s = [(2, 6), (11, 6), (2, 9), (11, 9), (2, 10), (11, 10), (8, 7), (11, 7), (5, 8), (11, 8)]
        alpha = [(0, 0), (0, 3), (3, 7), (11, 11), (2, 9), (7, 4)]
        beta = [(0, 0), (1, 1), (4, 0), (5, 4), (4, 3), (3, 5)]
    else:
        raise ValueError("code must be 'gross' or 'two-gross'.")
    return l, m, p, q, r, s, alpha, beta


def _align_z_rows_to_pairing(X_ops: np.ndarray, Z_ops: np.ndarray) -> np.ndarray:
    """Transform Z rows so X_ops @ Z_ops.T becomes identity."""
    P = gf2_matmul(X_ops, Z_ops.T)
    k = P.shape[0]
    if np.array_equal(P, np.eye(k, dtype=np.uint8)):
        return Z_ops
    if gf2_rank(P) != k:
        raise RuntimeError("Cannot align Z rows: pairing is singular.")
    # If P is a permutation, just reorder rows and preserve paper weights.
    if np.all(P.sum(axis=0) == 1) and np.all(P.sum(axis=1) == 1):
        perm = [int(np.flatnonzero(P[i])[0]) for i in range(k)]
        return Z_ops[perm].copy()
    return gf2_matmul(gf2_inverse_square(P).T, Z_ops)


def tour_de_gross_basis(code: Literal["gross", "two-gross"], *, reduce_weight: bool = False) -> LogicalBasis:
    """Paper-aligned logical basis for gross or two-gross.

    Uses A=1+y+x^3 y^-1, B=1+x+x^-1 y^-3 and the Appendix A.1 p,q,r,s,
    alpha,beta data.  Labels q1..q12 correspond to the paper's logical qubits.
    """
    l, m, p, q, r, s, alpha, beta = _tour_polys(code)
    A_norm = normalize_terms(_TOUR_A, l, m)
    B_norm = normalize_terms(_TOUR_B, l, m)
    qT, pT = transpose_terms(q, l, m), transpose_terms(p, l, m)
    sT, rT = transpose_terms(s, l, m), transpose_terms(r, l, m)
    mu = nu = (1, 1)
    X_rows: List[np.ndarray] = []
    Z_rows: List[np.ndarray] = []

    # Eq. (34), first block 1..6.
    for a, b in zip(alpha, beta):
        X_rows.append(lr_vec(shift_terms(p, a, l, m), shift_terms(q, a, l, m), l, m))
        bnu = (b[0] + nu[0], b[1] + nu[1])
        Z_rows.append(lr_vec(shift_terms(sT, bnu, l, m), shift_terms(rT, bnu, l, m), l, m))
    # Eq. (34), second block 7..12.
    for a, b in zip(alpha, beta):
        bT = (-b[0], -b[1])
        aTmu = (-a[0] + mu[0], -a[1] + mu[1])
        X_rows.append(lr_vec(shift_terms(r, bT, l, m), shift_terms(s, bT, l, m), l, m))
        Z_rows.append(lr_vec(shift_terms(qT, aTmu, l, m), shift_terms(pT, aTmu, l, m), l, m))

    X_ops = np.vstack(X_rows).astype(np.uint8)
    Z_ops = _align_z_rows_to_pairing(X_ops, np.vstack(Z_rows).astype(np.uint8))
    _, _, HX, HZ, _ = bb.build_blocks(l, m, A_norm, B_norm)
    if reduce_weight:
        # Usually keep this off: the paper basis is chosen for LPU geometry.
        Z_ops, X_ops = optimize_logicals_by_stabilizers(Z_ops, X_ops, HX, HZ)
    basis = LogicalBasis(
        l=l, m=m, A_exp=A_norm, B_exp=B_norm, Z_ops=Z_ops, X_ops=X_ops, HX=HX, HZ=HZ,
        labels=[f"q{i}" for i in range(1, 13)],
        source=f"Tour-de-gross Appendix A.1 {code} basis",
        notes="Named q1..q12 architecture basis; stabilizer reduction is disabled by default to preserve LPU/connectivity structure.",
    )
    if not basis.validate()["valid"]:
        raise RuntimeError(f"Tour-de-gross {code} basis failed validation: {basis.validate()}")
    return basis


def is_tour_gross_code(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> bool:
    """Return whether the input is the 12-by-6 Tour-de-gross BB code."""
    return (int(l), int(m)) == (12, 6) and set(normalize_terms(A_exp, l, m)) == set(normalize_terms(_TOUR_A, l, m)) and set(normalize_terms(B_exp, l, m)) == set(normalize_terms(_TOUR_B, l, m))


def is_tour_two_gross_code(l: int, m: int, A_exp: Iterable[Sequence[int]], B_exp: Iterable[Sequence[int]]) -> bool:
    """Return whether the input is the 12-by-12 two-gross BB code."""
    return (int(l), int(m)) == (12, 12) and set(normalize_terms(A_exp, l, m)) == set(normalize_terms(_TOUR_A, l, m)) and set(normalize_terms(B_exp, l, m)) == set(normalize_terms(_TOUR_B, l, m))


def best_logical_basis(
    l: int,
    m: int,
    A_exp: Iterable[Sequence[int]],
    B_exp: Iterable[Sequence[int]],
    *,
    prefer_tour_basis: bool = True,
    reduce_weight: bool = True,
) -> LogicalBasis:
    """Choose paper basis for gross/two-gross when applicable, else canonical."""
    A_norm = normalize_terms(A_exp, l, m)
    B_norm = normalize_terms(B_exp, l, m)
    if prefer_tour_basis:
        if is_tour_gross_code(l, m, A_norm, B_norm):
            return tour_de_gross_basis("gross", reduce_weight=False)
        if is_tour_two_gross_code(l, m, A_norm, B_norm):
            return tour_de_gross_basis("two-gross", reduce_weight=False)
    return canonical_logical_basis(l, m, A_norm, B_norm, reduce_weight=reduce_weight)


# ---------------------------------------------------------------------------
# Built-in cases, CLI, self-test
# ---------------------------------------------------------------------------


def builtin_case(name: str) -> Tuple[int, int, List[Exponent], List[Exponent]]:
    """Look up lattice dimensions and BB polynomial exponents by case name."""
    name = name.lower().replace("_", "-")
    if name in {"gross", "tour-gross"}:
        return 12, 6, list(_TOUR_A), list(_TOUR_B)
    if name in {"two-gross", "twogross", "tour-two-gross"}:
        return 12, 12, list(_TOUR_A), list(_TOUR_B)
    if name in {"144-standard", "bb144"}:
        return 12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]
    if name in {"288-standard", "bb288"}:
        return 12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)]
    if name in {"756", "bb756"}:
        return 21, 18, [(3, 0), (0, 10), (0, 17)], [(0, 5), (3, 0), (19, 0)]
    raise ValueError(f"Unknown built-in case: {name}")


def selftest() -> None:
    """Run built-in regression checks for dimensions and basis validity."""
    cases = [
        ("72", 6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
        ("90", 15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)], 8),
        ("108", 9, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 8),
        ("144", 12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
        ("288", 12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)], 12),
        ("360", 30, 6, [(9, 0), (0, 1), (0, 2)], [(0, 3), (25, 0), (26, 0)], 12),
        ("756", 21, 18, [(3, 0), (0, 10), (0, 17)], [(0, 5), (3, 0), (19, 0)], 16),
    ]
    for name, l, m, A, B, k in cases:
        basis = canonical_logical_basis(l, m, A, B, reduce_weight=True)
        rep = basis.validate()
        assert rep["valid"], name
        assert basis.k == k, name
        assert bb.k_bb_code_final(l, m, A, B) == k, name
        print(f"canonical {name:>3}: k={basis.k}, max_wX={max(rep['X_weights'])}, max_wZ={max(rep['Z_weights'])}")
    for code, wt in [("gross", 12), ("two-gross", 20)]:
        basis = tour_de_gross_basis(code)
        rep = basis.validate()
        assert rep["valid"] and all(w == wt for w in rep["X_weights"]) and all(w == wt for w in rep["Z_weights"]), code
        print(f"paper {code}: k={basis.k}, wt={wt}")
    diag = diagnose_legacy_logicals(12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)])
    assert diag["old_z_commutes_with_x_checks"] and diag["old_x_commutes_with_z_checks"]
    print("legacy 288-standard diagnostic:", {k: v for k, v in diag.items() if k != "old_pairing"})


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point for building and exporting one logical basis."""
    parser = argparse.ArgumentParser(description="Build paired logical operators for BB CSS codes.")
    parser.add_argument("--case", default="gross", help="gross, two-gross, 144-standard, 288-standard, 756")
    parser.add_argument("--json", default="", help="Write JSON report.")
    parser.add_argument("--markdown", default="", help="Write Markdown report.")
    parser.add_argument("--selftest", action="store_true", help="Run built-in tests.")
    parser.add_argument("--no-reduce", action="store_true", help="Disable stabilizer weight cleanup for generic canonical bases.")
    args = parser.parse_args(argv)
    if args.selftest:
        selftest()
        return 0
    l, m, A, B = builtin_case(args.case)
    basis = best_logical_basis(l, m, A, B, reduce_weight=not args.no_reduce)
    rep = basis.validate()
    print(f"case={args.case} l={l} m={m} n={basis.n} k={basis.k} source={basis.source}")
    print("valid=", rep["valid"], "canonical_pairing=", rep["canonical_pairing"])
    print("X_weights=", rep["X_weights"])
    print("Z_weights=", rep["Z_weights"])
    if args.json:
        basis.write_json(args.json, include_entries=True)
        print("wrote", args.json)
    if args.markdown:
        basis.write_markdown(args.markdown)
        print("wrote", args.markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
