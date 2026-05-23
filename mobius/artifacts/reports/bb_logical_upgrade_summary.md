# Logical-basis upgrade summary

## gross paper
- new basis source: Tour-de-gross Appendix A.1 gross basis
- n=144, k=12
- validation: valid=True, canonical_pairing=True, X quotient rank=12, Z quotient rank=12
- new X weights: [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
- new Z weights: [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
- legacy rows: Z=6, X=6, expected k=12, legacy pairing rank=0, legacy Z quotient increment=4, legacy X quotient increment=4

## two-gross paper
- new basis source: Tour-de-gross Appendix A.1 two-gross basis
- n=288, k=12
- validation: valid=True, canonical_pairing=True, X quotient rank=12, Z quotient rank=12
- new X weights: [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
- new Z weights: [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
- legacy rows: Z=6, X=6, expected k=12, legacy pairing rank=0, legacy Z quotient increment=0, legacy X quotient increment=0

## 144 standard
- new basis source: canonical quotient basis with IFFT seeds
- n=144, k=12
- validation: valid=True, canonical_pairing=True, X quotient rank=12, Z quotient rank=12
- new X weights: [48, 46, 46, 44, 40, 40, 48, 40, 32, 42, 44, 40]
- new Z weights: [36, 36, 36, 36, 24, 26, 26, 24, 24, 26, 24, 26]
- legacy rows: Z=6, X=6, expected k=12, legacy pairing rank=0, legacy Z quotient increment=4, legacy X quotient increment=4

## 288 standard
- new basis source: canonical quotient basis with IFFT seeds
- n=288, k=12
- validation: valid=True, canonical_pairing=True, X quotient rank=12, Z quotient rank=12
- new X weights: [76, 74, 66, 84, 64, 78, 74, 72, 76, 84, 78, 84]
- new Z weights: [50, 38, 38, 38, 38, 48, 48, 50, 38, 38, 44, 50]
- legacy rows: Z=6, X=6, expected k=12, legacy pairing rank=0, legacy Z quotient increment=0, legacy X quotient increment=0

## 756 standard
- new basis source: canonical quotient basis with IFFT seeds
- n=756, k=16
- validation: valid=True, canonical_pairing=True, X quotient rank=16, Z quotient rank=16
- new X weights: [182, 210, 180, 206, 192, 204, 174, 196, 188, 188, 192, 188, 188, 192, 188, 188]
- new Z weights: [188, 188, 192, 188, 188, 192, 188, 188, 168, 188, 198, 194, 194, 190, 168, 168]
- legacy rows: Z=8, X=8, expected k=16, legacy pairing rank=0, legacy Z quotient increment=8, legacy X quotient increment=8
