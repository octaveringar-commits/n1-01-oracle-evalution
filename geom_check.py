"""N1 geometry-only check (no channel draws): epochs without any admissible relay, per seed."""
import sys, csv, importlib.util
import numpy as np
n1dir, runs = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("n1_oracle", n1dir + "/n1_oracle.py")
n1 = importlib.util.module_from_spec(spec); sys.modules["n1_oracle"] = n1; spec.loader.exec_module(n1)
cfg = n1.N1Config()
out = {}
for r in csv.DictReader(open(runs)):
    if r["config"] == "0" and r["policy"] == "statique":
        out[int(r["seed"])] = round(float(r["outage"]) * cfg.episode_epochs)
print("seed  no_admissible_relay  archived_outage_epochs  match")
tot = 0
for s in range(20):
    ch = n1.N1Channel(cfg, s)
    none = sum(1 for t in range(cfg.episode_epochs) if not ch.geometry(t * cfg.dt_s)["vis"].any())
    tot += none
    print(f"{s:4d}  {none:19d}  {out[s]:22d}  {none == out[s]}")
print("total", tot, "archived", sum(out.values()))
