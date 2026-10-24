# N1-01 V1.3.3 — Clairvoyant Oracle Evaluation of Joint Relay–MCS Adaptation

## Description

Preregistered exploratory evaluation of joint relay selection and LDPC coding-rate adaptation in a Ground–LEO–GEO free-space optical architecture (NFIRE / EDRS terminal envelopes).

The campaign compares six policies on identical channel realizations (4 configurations × 20 seeds × 1000 epochs):

| Policy (code name) | Relay | MCS |
|---|---|---|
| `causale` (causal model-based) | argmax of predicted goodput Q_t | argmax of Q_t |
| `clairvoyante` (oracle, upper bound) | best realized action | best realized action |
| `r_fixe` (fixed relay) | highest-elevation relay at current geometry | causal (argmax_m Q_t) |
| `m_fixe` (fixed MCS) | argmax_R Q_t | constant, chosen ex ante |
| `statique` (static reference) | highest-elevation relay at current geometry | constant, chosen ex ante (argmax of Σ_t Q_t) |
| `directe` (direct OGS–GEO, no relay) | — | argmax_m of direct-link Q_t |

## Result (preregistered outcome: Issue C)

- All five relayed policies — including the clairvoyant oracle — yield exactly the same goodput and outage on all 80 trajectories (per-trajectory difference = 0 in floating point): J = 1.516239, P_out = 0.06405.
- Preregistered contrasts ΔJ_R, ΔJ_M and clairvoyant-vs-R: 0 in all 4 configurations (not significant).
- Preregistered contrast ΔJ_tot (relayed vs. direct): −6.4 % in all 4 configurations (95 % BCa interval entirely negative). The direct OGS–GEO link (J = 1.620, 0 % outage) outperforms every relayed policy.
- Conclusion: under the evaluated terminal envelopes, joint relay–MCS adaptation has no decision value.

## Files

| File | Content |
|---|---|
| `n1_oracle.py` | Evaluation code (frozen) |
| `test_n1_oracle.py` | Unit tests (formulas, invariants, governance) |
| `fso_env.py` | V19 revision 2 environment, imported unmodified (hash checked at import) |
| `runs.csv` | Per-run results (4 configurations × 20 seeds × 6 policies) |
| `contrasts.csv` | Preregistered contrasts, BCa intervals, exact sign-flip p-values |
| `summary.json` | Campaign outcome, hashes, significance flags |
| `SHA256SUMS.txt` | SHA-256 of runs.csv, contrasts.csv, summary.json |

## SHA-256

| File | SHA-256 |
|---|---|
| n1_oracle.py | `248b3dcf387043396881d2dc96612416577de0ce8a6386cd325ba6f2f8e32e3f` |
| fso_env.py | `a8e6dcbb77b0a21b40530cf285b455b667dc8c3e9f6ba685de69376135c3f851` |
| test_n1_oracle.py | `9f7ae5962bdb0b1f8494bfd1ebf1de00c68c16569c8a3c380af31a3c6f87fbd0` |
| runs.csv | `1ae486e925ad99617a2e030d10e6d629ce44159d5b59ecfe47290f8d5616ce29` |
| contrasts.csv | `944e28744213adf6f4c59cc7c8b4d78a25633954e32a155dd37a63bce6007617` |
| summary.json | `75c270829ad241d151be2d5e21e3c322228daef9e443a39778a24334a2089b59` |
| ETSI EN 302 307-1 V1.4.1 PDF (not distributed) | `19077f80420eb17519fcf6f83e2847d7a877940b134c5306814b8c6f60d8475d` |

## Reproduction

The frozen code expects `fso_env.py` in a sibling folder named `v19`. Recreate this layout without modifying any file:

```bash
mkdir -p n1 v19
cp n1_oracle.py test_n1_oracle.py n1/
cp fso_env.py v19/
cd n1 && python test_n1_oracle.py      # unit tests only, no grid run
sha256sum -c ../SHA256SUMS.txt         # run from the folder containing the result files
```

Requirements: Python 3, NumPy, SciPy. Running the full grid requires the explicit authorization token defined in the source.

## Scope

Results apply only to the evaluated terminal envelopes (NFIRE 5.625 GBd, EDRS 1.8 GBd), clear-sky assumption (τ = 0.10), untracked uplink, and decoding on the epoch average of 64 subslots. They do not claim a universal absence of adaptation value.

## License

CC0-1.0
