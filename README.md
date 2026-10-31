# N1-01 V1.3.3 — Clairvoyant Oracle Evaluation of Joint Relay–MCS Adaptation

Preregistered evaluation of joint relay selection and coding-rate adaptation in a Ground–LEO–GEO decode-and-forward free-space optical link model.

## Contents

| File | Description |
|---|---|
| `fso_env.py` | Simulation environment (V19 rev. 2) |
| `n1_oracle.py` | N1-01 V1.3.3 evaluation code |
| `test_n1_oracle.py` | Unit checks |
| `runs.csv` | Per-trajectory results |
| `contrasts.csv` | Per-contrast statistics |
| `summary.json` | Outcome, code and standard checksums |
| `SHA256SUMS.txt` | SHA-256 checksums of the campaign outputs |
| `LICENSE` | CC0-1.0 |

All files are at the repository root.

## Model (summary)

- Architecture: ground station → LEO relay → GEO, decode-and-forward; direct ground → GEO link as a separate comparator.
- Modulation: coherent homodyne BPSK at 1064 nm; nominal symbol rates 5.625 GBd (hop 1) and 1.8 GBd (hop 2).
- Coding: eleven DVB-S2 LDPC rates, one rate common to both hops; goodput ceiling 0.9 × 1.8 = 1.620 Gbit/s.
- Decoding thresholds: QPSK thresholds of ETSI EN 302 307-1 V1.4.1, Table 13, minus 3.01 dB (assumption).
- Receiver: 7.85 dB implementation loss (48 photons/bit at BER 1e-8); Es/N0 = 2·n_s/L.
- Uplink channel: Andrews et al. (2006), untracked case; long-term beam W_LT² = W² + ⟨r_c²⟩; log-normal intensity; 2 µrad mechanical jitter.
- Decision epoch: decoding on the mean of 64 sub-slots.
- Scope: clear sky, single ground station, uplink zenith angle ≤ 60°.

## Policies

| Policy | Relay | Coding rate |
|---|---|---|
| Clairvoyant | Best for the realized epoch channel | Best for the realized epoch channel (upper bound) |
| Causal | Joint best (R, M) on Q_t = E[G_t \| geometry] | Joint with relay |
| R-fixed | Highest-elevation visible relay | Best causal rate |
| M-fixed | Best causal relay | Constant per trajectory, chosen ex ante on ΣQ_t |
| Static | Highest-elevation visible relay | Constant per trajectory, chosen ex ante on ΣQ_t |
| Direct | Ground → GEO, no relay | Same rate set |

## Campaign

- Protocol: N1-01 V1.3.3, preregistered before execution.
- Configurations: W₀ ∈ {2.3, 3.25} cm × C_r ∈ {π, 3.86}.
- Design: 4 configurations × 20 seeds × 1,000 epochs (80 trajectories, 80,000 epochs).
- Statistics: exact sign-flip test, BCa intervals, Holm correction over the 4 configurations, ≥ 5 % relative threshold.
- Execution date: 21 October 2026.

## Outcome

- Preregistered outcome: **C** — no significant joint-adaptation gain under the preregistered criteria.
- All five relayed policies: mean goodput 1.516239 Gbit/s, outage fraction 0.06405.
- Direct ground-to-GEO comparator: 1.620 Gbit/s, no outage.
- The total contrast (−6.4 %) is a descriptive estimate; no reverse-direction test was preregistered.
- No general claim is made about adaptation in optical satellite networks.

## Integrity

All files are byte-identical to those of the 21 October 2026 campaign.

    sha256sum -c SHA256SUMS.txt
    sha256sum n1_oracle.py test_n1_oracle.py fso_env.py

Expected SHA-256:

    248b3dcf387043396881d2dc96612416577de0ce8a6386cd325ba6f2f8e32e3f  n1_oracle.py
    9f7ae5962bdb0b1f8494bfd1ebf1de00c68c16569c8a3c380af31a3c6f87fbd0  test_n1_oracle.py
    a8e6dcbb77b0a21b40530cf285b455b667dc8c3e9f6ba685de69376135c3f851  fso_env.py
    1ae486e925ad99617a2e030d10e6d629ce44159d5b59ecfe47290f8d5616ce29  runs.csv
    944e28744213adf6f4c59cc7c8b4d78a25633954e32a155dd37a63bce6007617  contrasts.csv
    75c270829ad241d151be2d5e21e3c322228daef9e443a39778a24334a2089b59  summary.json

## Execution lock

`run_grid()` requires the V1.3.3 authorization token and refuses a second run in the same directory. These safeguards were part of the preregistered protocol and are retained unchanged.

## Not included

- ETSI EN 302 307-1 V1.4.1 and Andrews et al. (2006), *Opt. Eng.* 45(7) 076001: not redistributed for copyright reasons; their SHA-256 checksums are recorded in `summary.json`.
- The post hoc link-margin diagnostic, which lies outside the preregistered protocol, is not part of this archive.

## Archive

Zenodo v1.0.0 (incomplete): https://doi.org/10.5281/zenodo.23612900. This release (v1.0.1) supersedes it.



## Post hoc diagnostics

The three scripts below were added in release v1.0.2. They lie **outside the prespecified protocol**: they do not modify the archived code, the campaign outputs or the contrasts, and they draw no new campaign result. They support Sections 3.4 and 3.5 of the manuscript.

| Script | Purpose | SHA-256 |
|---|---|---|
| `geom_check.py` | Counts, for each seed, the epochs without any admissible relay (geometry only, no channel draw) and compares them with the archived outage counts in `runs.csv` | `5a6c386096773141ed0833ff5b6270c19829d6b7a7c84f9e6b087aebd326c00a` |
| `seed15.py` | For seeds 2, 12 and 15: minimum zenith angle of the closest relay over the 100 s window, and number of epochs in which a relay lies within 60° but line of sight is blocked | `b4fd9398d9af0297d9d1352b48811d27edc1d02719658646088ef92fdaaea9b6` |
| `margins.py` | Link margins: 0.1 % lower quantile of the epoch-mean E_s/N_0 relative to the rate-9/10 threshold (seed 0, 2,000 epochs, four configurations) | `2cf492a94198c694546fafc35dc7c16222520d9a53ff317d8c33f0d78bd46bf8` |

Usage (`n1_oracle.py` expects `fso_env.py` in a sibling `v19/` folder):

```
mkdir -p run/n1 run/v19
cp n1_oracle.py run/n1/ && cp fso_env.py run/v19/
python geom_check.py run/n1 runs.csv
python seed15.py run/n1
python margins.py run/n1
```

Expected results: archived outages equal epochs without an admissible relay on every seed (158 for seed 2, 123 for seed 12, 1,000 for seed 15, 0 otherwise; 1,281 in total); uplink margins 30.7–40.6 dB over zenith angles 0–60°, inter-satellite margins 23.5–23.6 dB, direct-link margins 7.5–9.6 dB.

## License

CC0-1.0.
