"""
N1-01 post hoc link-margin diagnostic (outside the preregistered protocol).

Margin = 0.1 % lower quantile of the epoch-mean E_s/N0 (64 sub-slots, linear mean)
minus the BPSK threshold of the highest LDPC rate (9/10), in dB.
Uses the archived channel model of n1_oracle.py (V1.3.3) unchanged; seed 0, 2,000 epochs
per point; the four preregistered configurations (W0 x C_r).

  hop 1  (ground -> LEO, 550 km): zenith angles 0-60 deg, step 5 deg
  hop 2  (LEO -> GEO, vacuum):    range of LEO-GEO distances met by admissible relays (20 seeds)
  direct (ground -> GEO):         station-to-GEO geometry (zenith 14.3 deg)

Usage: python margins.py <dir containing n1_oracle.py; ../v19/fso_env.py next to it>
"""
import copy
import importlib.util
import sys

import numpy as np

spec = importlib.util.spec_from_file_location("n1_oracle", sys.argv[1] + "/n1_oracle.py")
n1 = importlib.util.module_from_spec(spec)
sys.modules["n1_oracle"] = n1
spec.loader.exec_module(n1)

Q, N_EP, SEED = 0.001, 2000, 0
thr_db = 10 * np.log10(n1.bpsk_thresholds_lin()[-1])
base = n1.N1Config()
K = base.n_subslots


def margin_db(esn0_lin):
    return 10 * np.log10(np.quantile(esn0_lin.mean(axis=1), Q)) - thr_db


# LEO-GEO distance range met by admissible relays (geometry only, 20 seeds)
d2 = []
for s in range(20):
    ch = n1.N1Channel(base, s)
    for t in range(base.episode_epochs):
        g = ch.geometry(t * base.dt_s)
        d2.extend(g["d2"][g["vis"]].tolist())
d2_min, d2_max = min(d2), max(d2)
zen_geo = n1.N1Channel(base, 0).geometry(0.0)["zen_geo"]

print(f"threshold r=9/10 (BPSK): {thr_db:.2f} dB | GEO zenith {np.degrees(zen_geo):.2f} deg | "
      f"LEO-GEO distance {d2_min/1e3:.0f}-{d2_max/1e3:.0f} km")
print("config  W0[cm]  C_r    hop1 margin by zenith 0..60 deg (dB)                                  hop2 (dB)       direct (dB)")
for k, params in enumerate(n1.GRID):
    cfg = copy.deepcopy(base)
    for key, val in params.items():
        setattr(cfg, key, val)
    cm = n1.ChannelModel(cfg)
    rng = np.random.default_rng(SEED)
    h1 = [margin_db(cm.uplink_esn0(np.radians(z), cfg.v19.h_leo_m, cfg.rs1_baud, cfg.d_rx_leo_m, (N_EP, K), rng))
          for z in range(0, 61, 5)]
    h2 = [margin_db(cm.isl_esn0(d, (N_EP, K), rng)) for d in (d2_min, d2_max)]
    hd = margin_db(cm.uplink_esn0(zen_geo, cfg.v19.h_geo_m, cfg.rs2_baud, cfg.d_rx_geo_m, (N_EP, K), rng))
    print(f"{k:6d}  {cfg.w0_gs_m*100:6.2f}  {cfg.c_r:4.2f}  " + " ".join(f"{m:5.1f}" for m in h1)
          + f"   {h2[0]:5.1f}-{h2[1]:5.1f}   {hd:5.1f}")
