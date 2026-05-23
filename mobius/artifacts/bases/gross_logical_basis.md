# Logical basis report: Tour-de-gross Appendix A.1 gross basis

- l=12, m=6, n=144, k=12
- A exponents: [(0, 0), (0, 1), (3, 5)]
- B exponents: [(0, 0), (1, 0), (11, 3)]
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
- X_weights: [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
- Z_weights: [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]

## Logical qubits
| logical | wt(X) | X support | wt(Z) | Z support |
|---|---:|---|---:|---|
| q1 | 12 | `L: x^4 + x^4y^2 + x^5 + x^5y^4 + x^6y + x^6y^5 | R: x^3 + x^3y + x^3y^2 + x^3y^5 + x^4 + x^4y^2` | 12 | `L: y + xy^2 + x^4y + x^5y^2 + x^5y^3 + x^9y^3 | R: 1 + xy + x^2y^3 + x^4 + x^5y + x^10y^3` |
| q2 | 12 | `L: x^7y + x^7y^5 + x^8y^3 + x^8y^5 + x^9 + x^9y^4 | R: x^6 + x^6y + x^6y^4 + x^6y^5 + x^7y + x^7y^5` | 12 | `L: x^2y^5 + x^3 + x^6y^5 + x^7 + x^7y + x^11y | R: y + x^2y^4 + x^3y^5 + x^4y + x^6y^4 + x^7y^5` |
| q3 | 12 | `L: x^3y + x^3y^5 + x^4y^3 + x^4y^5 + x^5 + x^5y^4 | R: x^2 + x^2y + x^2y^4 + x^2y^5 + x^3y + x^3y^5` | 12 | `L: xy^3 + x^2y^4 + x^5y^3 + x^6y^4 + x^6y^5 + x^10y^5 | R: xy^2 + x^2y^3 + x^3y^5 + x^5y^2 + x^6y^3 + x^11y^5` |
| q4 | 12 | `L: x^2y + x^2y^3 + x^3y + x^3y^5 + x^4 + x^4y^2 | R: x + xy + xy^2 + xy^3 + x^2y + x^2y^3` | 12 | `L: x^2 + x^3y + x^6 + x^7y + x^7y^2 + x^11y^2 | R: y^2 + x^2y^5 + x^3 + x^4y^2 + x^6y^5 + x^7` |
| q5 | 12 | `L: x^9 + x^9y^4 + x^10y^2 + x^10y^4 + x^11y^3 + x^11y^5 | R: x^8 + x^8y^3 + x^8y^4 + x^8y^5 + x^9 + x^9y^4` | 12 | `L: xy^2 + x^2y^3 + x^5y^2 + x^6y^3 + x^6y^4 + x^10y^4 | R: xy + x^2y^2 + x^3y^4 + x^5y + x^6y^2 + x^11y^4` |
| q6 | 12 | `L: x^8y^2 + x^8y^4 + x^9 + x^9y^2 + x^10y + x^10y^3 | R: x^7y + x^7y^2 + x^7y^3 + x^7y^4 + x^8y^2 + x^8y^4` | 12 | `L: y^4 + x^3y^2 + x^4y^3 + x^7y^2 + x^8y^3 + x^8y^4 | R: xy^4 + x^3y + x^4y^2 + x^5y^4 + x^7y + x^8y^2` |
| q7 | 12 | `L: 1 + xy + x^3y^4 + x^8 + x^9y + x^11y^4 | R: y^5 + x + x^4y^4 + x^8y^4 + x^8y^5 + x^9` | 12 | `L: x^9y + x^9y^5 + x^10 + x^10y + x^10y^2 + x^10y^5 | R: x^7 + x^7y^2 + x^8y + x^8y^3 + x^9y + x^9y^5` |
| q8 | 12 | `L: x + x^6y^2 + x^7y^3 + x^9 + x^10y^2 + x^11y^3 | R: x^2 + x^6 + x^6y + x^7y^2 + x^10y + x^11y^2` | 12 | `L: x^6 + x^6y^2 + x^7 + x^7y + x^7y^2 + x^7y^3 | R: x^4y + x^4y^3 + x^5y^2 + x^5y^4 + x^6 + x^6y^2` |
| q9 | 12 | `L: y^5 + x^2y^2 + x^7y^4 + x^8y^5 + x^10y^2 + x^11y^4 | R: y^4 + x^3y^2 + x^7y^2 + x^7y^3 + x^8y^4 + x^11y^3` | 12 | `L: x^10 + x^10y^2 + x^11 + x^11y + x^11y^2 + x^11y^3 | R: x^8y + x^8y^3 + x^9y^2 + x^9y^4 + x^10 + x^10y^2` |
| q10 | 12 | `L: xy^5 + x^6y + x^7y^2 + x^9y^5 + x^10y + x^11y^2 | R: x^2y^5 + x^6 + x^6y^5 + x^7y + x^10 + x^11y` | 12 | `L: 1 + y + y^4 + y^5 + x^11 + x^11y^4 | R: x^9y + x^9y^5 + x^10 + x^10y^2 + x^11 + x^11y^4` |
| q11 | 12 | `L: 1 + x^2y^3 + x^7y^5 + x^8 + x^10y^3 + x^11y^5 | R: y^5 + x^3y^3 + x^7y^3 + x^7y^4 + x^8y^5 + x^11y^4` | 12 | `L: x^4y + x^4y^3 + x^5y + x^5y^2 + x^5y^3 + x^5y^4 | R: x^2y^2 + x^2y^4 + x^3y^3 + x^3y^5 + x^4y + x^4y^3` |
| q12 | 12 | `L: y^3 + x^5y^5 + x^6 + x^8y^3 + x^9y^5 + x^10 | R: xy^3 + x^5y^3 + x^5y^4 + x^6y^5 + x^9y^4 + x^10y^5` | 12 | `L: x^5y^3 + x^5y^5 + x^6 + x^6y^3 + x^6y^4 + x^6y^5 | R: x^3 + x^3y^4 + x^4y + x^4y^5 + x^5y^3 + x^5y^5` |
