import sys, importlib.util
import numpy as np
spec = importlib.util.spec_from_file_location("n1_oracle", sys.argv[1] + "/n1_oracle.py")
n1 = importlib.util.module_from_spec(spec); sys.modules["n1_oracle"] = n1; spec.loader.exec_module(n1)
v19 = n1.v19; cfg = n1.N1Config()
for s in (2, 12, 15):
    ch = n1.N1Channel(cfg, s)
    zmin, blocked = [], 0
    for t in range(cfg.episode_epochs):
        g = ch.geometry(t * cfg.dt_s)
        z = np.degrees(g["zen"]); zmin.append(z.min())
        within = z <= cfg.zenith_max_deg
        if within.any() and not g["vis"].any(): blocked += 1
    print(s, "min zenith over 100 s: %.1f deg" % min(zmin), "| epochs with relay <=60deg but LOS blocked:", blocked)
