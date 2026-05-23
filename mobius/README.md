# Mobius BB Logical Tools

This directory contains standalone utilities and generated artifacts for
Mobius-style bivariate-bicycle (BB) qLDPC code experiments.

## Layout

- `bb_fourier_tools_fast.py`: finite-field, Fourier-rank, GF(2) linear algebra,
  and canonical logical-basis helpers.
- `bb_logical_tools.py`: legacy paired-logical helpers used by the original
  smoke tests.
- `bb_logical_tools_clean.py`: maintained command-line and library interface for
  paper-aligned and generic paired logical bases.
- `tests/`: executable regression/smoke tests.
- `artifacts/bases/`: generated logical-basis JSON and Markdown reports.
- `artifacts/matrices/`: generated `.npz` logical matrix exports.
- `artifacts/pairs/`: generated per-logical-qubit support JSON files.
- `artifacts/reports/`: generated summaries and captured test output.

## Common Commands

Run the cleaned module's self-test:

```bash
python mobius/bb_logical_tools_clean.py --selftest
```

Build and export one named basis:

```bash
python mobius/bb_logical_tools_clean.py \
  --case gross \
  --json mobius/artifacts/bases/gross_logical_basis.json \
  --markdown mobius/artifacts/bases/gross_logical_basis.md
```

Run the Mobius smoke tests directly:

```bash
python mobius/tests/test_bb_logical_tools.py
python mobius/tests/test_bb_logical_tools_clean.py
```
