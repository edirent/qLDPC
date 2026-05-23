"""Smoke tests for the legacy paired-logical helper module.

The cases below exercise the standard BB examples used while developing the
Mobius logical-operator scripts.  They intentionally verify both the algebraic
dimension ``k`` and the exported support format, because downstream notebooks
consume those JSON artifacts directly.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

MOBIUS_ROOT = Path(__file__).resolve().parents[1]
if str(MOBIUS_ROOT) not in sys.path:
    sys.path.insert(0, str(MOBIUS_ROOT))

from bb_logical_tools import (
    logical_ops_bb_paired,
    validate_css_logicals,
    logical_qubits_from_ops,
    qldpc_style_logical_matrix,
    export_logical_qubits_json,
    summarize_logicals,
    tour_de_gross_logicals,
)
from bb_fourier_tools_fast import k_bb_code_final

STANDARD_CASES = [
    ("[[72,12,6]]", 6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
    ("[[90,8,10]]", 15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)], 8),
    ("[[108,8,10]]", 9, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 8),
    ("[[144,12,12]]", 12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
    ("[[288,12,18]]", 12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)], 12),
    ("[[360,12,<=24]]", 30, 6, [(9, 0), (0, 1), (0, 2)], [(0, 3), (25, 0), (26, 0)], 12),
    ("[[756,16,<=34]]", 21, 18, [(3, 0), (0, 10), (0, 17)], [(0, 5), (3, 0), (19, 0)], 16),
]

PAIRS_DIR = MOBIUS_ROOT / "artifacts" / "pairs"
REPORTS_DIR = MOBIUS_ROOT / "artifacts" / "reports"


def check_case(name, l, m, A, B, expected_k):
    """Build one BB basis, validate it, and optionally refresh JSON artifacts."""
    Z, X, HX, HZ = logical_ops_bb_paired(l, m, A, B, reduce=True, reduction_trials=0, seed=123)
    report = validate_css_logicals(HX, HZ, X, Z)
    k = k_bb_code_final(l, m, A, B)
    assert report["valid"], report
    assert k == expected_k == X.shape[0] == Z.shape[0]
    assert qldpc_style_logical_matrix(X, Z).shape == (2 * k, 2 * HX.shape[1])
    logicals = logical_qubits_from_ops(X, Z, l, m)
    assert len(logicals) == k
    if name == "[[144,12,12]]":
        export_logical_qubits_json(logicals, str(PAIRS_DIR / "bb_logical_pairs_144_generic.json"))
    return {"case": name, "summary": summarize_logicals(X, Z), "valid": report["valid"]}


def main():
    """Run every smoke case and refresh the compact summary artifact."""
    PAIRS_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for case in STANDARD_CASES:
        rows.append(check_case(*case))

    for code in ["gross", "two-gross"]:
        l, m, A, B, Z, X, HX, HZ = tour_de_gross_logicals(code)
        report = validate_css_logicals(HX, HZ, X, Z)
        assert report["valid"], report
        rows.append({"case": f"tour-{code}", "summary": summarize_logicals(X, Z), "valid": report["valid"]})
        logicals = logical_qubits_from_ops(X, Z, l, m)
        export_logical_qubits_json(logicals, str(PAIRS_DIR / f"bb_logical_pairs_tour_{code}.json"))

    with open(REPORTS_DIR / "bb_logical_tools_test_summary.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)
    for row in rows:
        print(row["case"], row["summary"])
    print("all logical-tool tests passed")


def test_standard_and_tour_logical_tools():
    """Pytest entry point that mirrors direct script execution."""
    main()


if __name__ == "__main__":
    main()
