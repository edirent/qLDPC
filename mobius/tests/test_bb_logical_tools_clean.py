"""Regression checks for the cleaned Mobius logical-basis module.

The script focuses on the maintained ``bb_logical_tools_clean`` API: named
Tour-de-gross bases, generic canonical bases, and the legacy-kernel diagnostic.
It can be run directly or collected by pytest.
"""

from __future__ import annotations

from pathlib import Path
import sys

MOBIUS_ROOT = Path(__file__).resolve().parents[1]
if str(MOBIUS_ROOT) not in sys.path:
    sys.path.insert(0, str(MOBIUS_ROOT))

import bb_logical_tools_clean as lt


def test_clean_logical_tools_smoke():
    """Exercise the maintained API on paper-aligned and generic BB cases."""
    # Paper-aligned bases keep the Appendix A.1 weights and labels.
    for case in ["gross", "two-gross"]:
        l, m, A, B = lt.builtin_case(case)
        basis = lt.best_logical_basis(l, m, A, B)
        rep = basis.validate()
        assert rep["valid"], case
        assert basis.qldpc_style_logical_matrix().shape == (2 * basis.k, 2 * basis.n)
        print(
            case,
            "paper basis ok",
            "k=",
            basis.k,
            "Xmax=",
            max(rep["X_weights"]),
            "Zmax=",
            max(rep["Z_weights"]),
        )

    # Generic canonical bases cover standard BB conventions not covered by the
    # paper's named gross/two-gross construction.
    for case in ["144-standard", "288-standard", "756"]:
        l, m, A, B = lt.builtin_case(case)
        basis = lt.canonical_logical_basis(l, m, A, B, reduce_weight=False)
        rep = basis.validate()
        assert rep["valid"], case
        print(case, "canonical basis ok", "k=", basis.k)

    # The old IFFT-kernel representatives commute with checks, but this
    # diagnostic verifies that they are not a full paired basis for 288-standard.
    l, m, A, B = lt.builtin_case("288-standard")
    diag = lt.diagnose_legacy_logicals(l, m, A, B)
    assert diag["old_z_commutes_with_x_checks"] and diag["old_x_commutes_with_z_checks"]
    assert not diag["old_rows_equal_k"]
    print("legacy diagnostic ok", {k: v for k, v in diag.items() if k != "old_pairing"})


if __name__ == "__main__":
    test_clean_logical_tools_smoke()
