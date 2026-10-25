# N1-01 V1.3.3 — Clairvoyant Oracle Evaluation of Joint Relay–MCS Adaptation

Preregistered evaluation of joint relay selection and coding-rate adaptation in a Ground–LEO–GEO decode-and-forward free-space optical link model (coherent BPSK, 1064 nm, eleven DVB-S2 LDPC rates).

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

## Campaign

- Protocol: N1-01 V1.3.3, preregistered before execution.
- Design: 4 configurations × 20 seeds × 1,000 epochs (80 trajectories).
- Execution date: 21 October 2026.

## Outcome

- Preregistered outcome: **C** — no significant joint-adaptation gain under the preregistered criteria.
- All five relayed policies: mean goodput 1.516239 Gbit/s, outage fraction 0.06405.
- Direct ground-to-GEO comparator: 1.620 Gbit/s, no outage.
- The total contrast (−6.4 %) is a descriptive estimate; no reverse-direction test was preregistered.
- Scope: nominal rates, clear sky, single ground station, zenith angle ≤ 60°. No general claim about adaptation in optical satellite networks.

## Integrity

All files are byte-identical to those of the 21 October 2026 campaign. To check:

    sha256sum -c SHA256SUMS.txt
    sha256sum n1_oracle.py test_n1_oracle.py fso_env.py

Expected SHA-256 prefixes: `n1_oracle.py` 248b3dcf, `test_n1_oracle.py` 9f7ae596, `fso_env.py` a8e6dcbb.

## Execution lock

`run_grid()` requires the V1.3.3 authorization token and refuses a second run in the same directory. These safeguards were part of the preregistered protocol and are retained unchanged.

## Not included

ETSI EN 302 307-1 V1.4.1 and Andrews et al. (2006), *Opt. Eng.* 45(7) 076001, are not redistributed for copyright reasons. Their SHA-256 checksums are recorded in `summary.json`.

## Archive

Zenodo: https://doi.org/10.5281/zenodo.23612900 (v1.0.0, incomplete). Release v1.0.1 supersedes it.

## License

CC0-1.0.
