# BB logical operator tools report

## What changed

- Replaced the old kernel-only logical construction with a complete paired CSS logical basis constructor.
- Returns exactly k logical X/Z pairs and validates X_ops @ Z_ops.T = I over GF(2).
- Adds optional stabilizer-coset weight reduction. This is heuristic and does not certify distance.
- Adds qLDPC-style symplectic output with shape (2*k, 2*n).
- Adds support-coordinate export: first l*m qubits are L, next l*m are R, with coordinates (x,y) on the torus.
- Adds Tour-de-gross Appendix A.1 bases for gross and two-gross in the paper convention.

## Legacy check

The original IFFT/kernel-style representatives return only dim(ker(A) ∩ ker(B)) rows per Pauli type, i.e. k/2 for BB codes.  For example, [[144,12,12]] and [[288,12,18]] return 6 X rows and 6 Z rows, while a full logical basis needs 12 X/Z pairs.  The new constructor outputs k pairs and checks each pair directly.

## Test summary

| case | k | n | X weight range | Z weight range | total weight |
|---|---:|---:|---:|---:|---:|
| [[72,12,6]] | 12 | 72 | 6–14 | 6–22 | 316 |
| [[90,8,10]] | 8 | 90 | 10–18 | 10–24 | 236 |
| [[108,8,10]] | 8 | 108 | 12–24 | 12–26 | 258 |
| [[144,12,12]] | 12 | 144 | 12–14 | 12–24 | 382 |
| [[288,12,18]] | 12 | 288 | 18–28 | 56–80 | 1082 |
| [[360,12,<=24]] | 12 | 360 | 44–94 | 30–66 | 1456 |
| [[756,16,<=34]] | 16 | 756 | 172–198 | 52–190 | 4782 |
| tour-gross | 12 | 144 | 12–12 | 12–12 | 288 |
| tour-two-gross | 12 | 288 | 20–20 | 20–20 | 480 |

All rows passed validation: commute with checks, nontrivial modulo stabilizers, and canonical pairings.