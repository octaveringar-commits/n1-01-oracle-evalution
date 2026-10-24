
Dr. octave ringar <octaveringar@gmail.com>
08:51 (il y a 5 minutes)
À moi

# N1-01 V1.3.3 — Clairvoyant Oracle Evaluation of Joint Relay–MCS Adaptation

## Description

Preregistered exploratory evaluation of joint relay selection and LDPC coding-rate adaptation in a Ground–LEO–GEO free-space optical architecture (NFIRE / EDRS terminal envelopes).

The campaign compares six policies on identical channel realizations:
- causal model-based
- clairvoyant oracle (upper bound)
- fixed-relay
- fixed-MCS
- static (highest-elevation relay + fixed MCS)
- direct OGS–GEO

**Result (Issue C):** no policy improves upon the static policy. The clairvoyant oracle attains exactly the same goodput as the static policy on all 80 archived trajectories.

## Files

| File | Content |
|------|---------|
| `n1_oracle.py` | Complete evaluation code (V19 revision 2 imported, frozen) |
| `runs.csv` | Aggregated results (4 configurations × 20 seeds × 6 policies) |
| `summary.json` | Campaign outcome, hashes, significance flags |

## Reproducibility

- Code SHA-256: `248b3dcf387043396881d2dc96612416577de0ce8a6386cd325ba6f2f8e32e3f`
- V19 SHA-256: `a8e6dcbb77b0a21b40530cf285b455b667dc8c3e9f6ba685de69376135c3f851`
- ETSI Table 13 SHA-256: `19077f80420eb17519fcf6f83e2847d7a877940b134c5306814b8c6f60d8475d`

Execution of the full grid requires the explicit authorization token defined in the source.

## Scope

Results apply only to the evaluated terminal envelopes (NFIRE 5.625 GBd, EDRS 1.8 GBd), clear-sky assumption (τ = 0.10), and decoding on the epoch average. They do not claim a universal absence of adaptation value.

## License

CC0-1.0
