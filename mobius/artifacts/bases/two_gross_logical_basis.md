# Logical basis report: Tour-de-gross Appendix A.1 two-gross basis

- l=12, m=12, n=288, k=12
- A exponents: [(0, 0), (0, 1), (3, 11)]
- B exponents: [(0, 0), (1, 0), (11, 9)]
- notes: Named q1..q12 architecture basis; stabilizer reduction is disabled by default to preserve LPU/connectivity structure.

## Validation
- valid: True
- expected_k: 12
- num_x: 12
- num_z: 12
- z_commutes_with_x_checks: True
- x_commutes_with_z_checks: True
- canonical_pairing: True
- z_nonstabilizer_quotient_rank: 12
- x_nonstabilizer_quotient_rank: 12
- pairing_rank: 12
- X_weights: [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
- Z_weights: [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]

## Logical qubits
| logical | wt(X) | X support | wt(Z) | Z support |
|---|---:|---|---:|---|
| q1 | 20 | `L: x^2 + x^2y^2 + x^3y^3 + x^4y^3 + x^6y^7 + x^7y^4 + x^7y^11 + x^8y^2 + x^8y^6 + x^9y^2 | R: xy^9 + xy^10 + x^3y^2 + x^4y^8 + x^5y^3 + x^6y^7 + x^7y^3 + x^8y^4 + x^8y^5 + x^11y^3` | 20 | `L: x^2y^3 + x^2y^4 + x^2y^5 + x^2y^6 + x^2y^7 + x^5y^6 + x^8y^5 + x^11y^3 + x^11y^4 + x^11y^7 | R: y^5 + y^8 + xy^8 + x^2y^11 + x^3y^2 + x^5y^5 + x^7y^8 + x^8y^8 + x^9y^11 + x^11y^2` |
| q2 | 20 | `L: x^2y^3 + x^2y^5 + x^3y^6 + x^4y^6 + x^6y^10 + x^7y^2 + x^7y^7 + x^8y^5 + x^8y^9 + x^9y^5 | R: x + xy + x^3y^5 + x^4y^11 + x^5y^6 + x^6y^10 + x^7y^6 + x^8y^7 + x^8y^8 + x^11y^6` | 20 | `L: y^4 + y^5 + y^8 + x^3y^4 + x^3y^5 + x^3y^6 + x^3y^7 + x^3y^8 + x^6y^7 + x^9y^6 | R: y^3 + xy^6 + xy^9 + x^2y^9 + x^3 + x^4y^3 + x^6y^6 + x^8y^9 + x^9y^9 + x^10` |
| q3 | 20 | `L: y^9 + x^5y^7 + x^5y^9 + x^6y^10 + x^7y^10 + x^9y^2 + x^10y^6 + x^10y^11 + x^11y + x^11y^9 | R: x^2y^10 + x^4y^4 + x^4y^5 + x^6y^9 + x^7y^3 + x^8y^10 + x^9y^2 + x^10y^10 + x^11 + x^11y^11` | 20 | `L: y^5 + x^3y^3 + x^3y^4 + x^3y^7 + x^6y^3 + x^6y^4 + x^6y^5 + x^6y^6 + x^6y^7 + x^9y^6 | R: y^8 + xy^11 + x^3y^2 + x^4y^5 + x^4y^8 + x^5y^8 + x^6y^11 + x^7y^2 + x^9y^5 + x^11y^8` |
| q4 | 20 | `L: xy + xy^11 + x^2y^2 + x^3y^2 + x^5y^6 + x^6y^3 + x^6y^10 + x^7y + x^7y^5 + x^8y | R: y^8 + y^9 + x^2y + x^3y^7 + x^4y^2 + x^5y^6 + x^6y^2 + x^7y^3 + x^7y^4 + x^10y^2` | 20 | `L: xy^9 + x^4y^7 + x^4y^8 + x^4y^11 + x^7y^7 + x^7y^8 + x^7y^9 + x^7y^10 + x^7y^11 + x^10y^10 | R: 1 + x + x^2y^3 + x^4y^6 + x^5 + x^5y^9 + x^6 + x^7y^3 + x^8y^6 + x^10y^9` |
| q5 | 20 | `L: x^4y^9 + x^4y^11 + x^5 + x^6 + x^8y^4 + x^9y + x^9y^8 + x^10y^3 + x^10y^11 + x^11y^11 | R: x + x^3y^6 + x^3y^7 + x^5y^11 + x^6y^5 + x^7 + x^8y^4 + x^9 + x^10y + x^10y^2` | 20 | `L: y^8 + x^3y^6 + x^3y^7 + x^3y^10 + x^6y^6 + x^6y^7 + x^6y^8 + x^6y^9 + x^6y^10 + x^9y^9 | R: y^11 + xy^2 + x^3y^5 + x^4y^8 + x^4y^11 + x^5y^11 + x^6y^2 + x^7y^5 + x^9y^8 + x^11y^11` |
| q6 | 20 | `L: xy^11 + x^2y^3 + x^2y^8 + x^3y^6 + x^3y^10 + x^4y^6 + x^9y^4 + x^9y^6 + x^10y^7 + x^11y^7 | R: y^7 + xy^11 + x^2y^7 + x^3y^8 + x^3y^9 + x^6y^7 + x^8y + x^8y^2 + x^10y^6 + x^11` | 20 | `L: x^2 + x^2y^8 + x^2y^9 + x^5 + x^5y^8 + x^5y^9 + x^5y^10 + x^5y^11 + x^8y^11 + x^11y^10 | R: y^4 + x^2y^7 + x^3y + x^3y^10 + x^4y + x^5y^4 + x^6y^7 + x^8y^10 + x^10y + x^11y` |
| q7 | 20 | `L: y^5 + xy^5 + xy^8 + x^2y^11 + x^4y^2 + x^5y^5 + x^6y^5 + x^8y^8 + x^10y^11 + x^11y^2 | R: x^2y^6 + x^2y^9 + x^2y^10 + x^5y^8 + x^8y^7 + x^11y^6 + x^11y^7 + x^11y^8 + x^11y^9 + x^11y^10` | 20 | `L: y^3 + y^4 + x^2y^10 + x^5y^8 + x^5y^9 + x^6y^10 + x^7y^6 + x^8y^10 + x^9y^5 + x^10y^11 | R: x^4y^11 + x^5y^7 + x^5y^11 + x^6y^2 + x^6y^9 + x^7y^6 + x^9y^10 + x^10y^10 + x^11y + x^11y^11` |
| q8 | 20 | `L: y^4 + y^7 + xy^10 + x^3y + x^4y^4 + x^5y^4 + x^7y^7 + x^9y^10 + x^10y + x^11y^4 | R: xy^5 + xy^8 + xy^9 + x^4y^7 + x^7y^6 + x^10y^5 + x^10y^6 + x^10y^7 + x^10y^8 + x^10y^9` | 20 | `L: 1 + y + x^2y^7 + x^5y^5 + x^5y^6 + x^6y^7 + x^7y^3 + x^8y^7 + x^9y^2 + x^10y^8 | R: x^4y^8 + x^5y^4 + x^5y^8 + x^6y^6 + x^6y^11 + x^7y^3 + x^9y^7 + x^10y^7 + x^11y^8 + x^11y^10` |
| q9 | 20 | `L: y^2 + xy^5 + x^2y^5 + x^4y^8 + x^6y^11 + x^7y^2 + x^8y^5 + x^9y^5 + x^9y^8 + x^10y^11 | R: xy^8 + x^4y^7 + x^7y^6 + x^7y^7 + x^7y^8 + x^7y^9 + x^7y^10 + x^10y^6 + x^10y^9 + x^10y^10` | 20 | `L: x^2y + x^2y^2 + x^3y^3 + x^4y^11 + x^5y^3 + x^6y^10 + x^7y^4 + x^9y^8 + x^9y^9 + x^11y^3 | R: xy^4 + x^2 + x^2y^4 + x^3y^2 + x^3y^7 + x^4y^11 + x^6y^3 + x^7y^3 + x^8y^4 + x^8y^6` |
| q10 | 20 | `L: y + xy + x^3y^4 + x^5y^7 + x^6y^10 + x^7y + x^8y + x^8y^4 + x^9y^7 + x^11y^10 | R: y^4 + x^3y^3 + x^6y^2 + x^6y^3 + x^6y^4 + x^6y^5 + x^6y^6 + x^9y^2 + x^9y^5 + x^9y^6` | 20 | `L: xy^4 + xy^5 + x^3y^11 + x^6y^9 + x^6y^10 + x^7y^11 + x^8y^7 + x^9y^11 + x^10y^6 + x^11 | R: 1 + y^2 + x^5 + x^6 + x^6y^8 + x^7y^3 + x^7y^10 + x^8y^7 + x^10y^11 + x^11y^11` |
| q11 | 20 | `L: y^11 + xy^2 + x^2y^2 + x^4y^5 + x^6y^8 + x^7y^11 + x^8y^2 + x^9y^2 + x^9y^5 + x^10y^8 | R: xy^5 + x^4y^4 + x^7y^3 + x^7y^4 + x^7y^5 + x^7y^6 + x^7y^7 + x^10y^3 + x^10y^6 + x^10y^7` | 20 | `L: y + x^3 + x^3y^11 + x^4y + x^5y^9 + x^6y + x^7y^8 + x^8y^2 + x^10y^6 + x^10y^7 | R: x^2y^2 + x^3y^2 + x^3y^10 + x^4 + x^4y^5 + x^5y^9 + x^7y + x^8y + x^9y^2 + x^9y^4` |
| q12 | 20 | `L: xy^9 + x^2 + x^3 + x^5y^3 + x^7y^6 + x^8y^9 + x^9 + x^10 + x^10y^3 + x^11y^6 | R: x^2y^3 + x^5y^2 + x^8y + x^8y^2 + x^8y^3 + x^8y^4 + x^8y^5 + x^11y + x^11y^4 + x^11y^5` | 20 | `L: y^2 + xy^6 + x^2y + x^3y^7 + x^5 + x^5y^11 + x^7y^6 + x^10y^4 + x^10y^5 + x^11y^6 | R: y^2 + x^2y^6 + x^3y^6 + x^4y^7 + x^4y^9 + x^9y^7 + x^10y^3 + x^10y^7 + x^11y^5 + x^11y^10` |
