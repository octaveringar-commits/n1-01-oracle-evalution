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

## License

CC0-1.0.
