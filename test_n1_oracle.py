"""
test_n1_oracle.py — Tests unitaires de n1_oracle.py (N1-01 V1.3.3).
Ne lance PAS la grille : seuls des contrôles de formules, d'invariants et de
gouvernance. Aucune valeur de J n'est imprimée.
"""
import copy

import numpy as np

import n1_oracle as n1

PASSED = []


def check(name, cond, detail=""):
    assert cond, f"ÉCHEC {name} {detail}"
    PASSED.append(name)
    print(f"  OK  {name} {detail}")


cfg = n1.N1Config()

print("A — Canal montant certifié (Andrews 2006)")
# A1 : reproduction de la Fig. 6 (1,6 µm, 300 km, HV5/7, zénith)
c6 = copy.deepcopy(cfg); c6.lambda_m = 1.6e-6
ref = {0.005: 0.052, 0.01: 0.045, 0.03: 0.027, 0.10: 0.008}
for w0, val in ref.items():
    c6.w0_gs_m = w0
    st = n1.uplink_stats(0.0, 300e3, c6)
    check(f"A1 Fig. 6 σ²_B, W0={w0*100:g} cm", abs(st["sigma_b2"] / val - 1) < 0.10,
          f"{st['sigma_b2']:.4f} vs ≈{val}")
check("A2 Fig. 6 r0 = 20 cm à 1,6 µm", abs(n1.uplink_stats(0.0, 300e3, c6)["r0"] - 0.20) < 0.005)
# A3 : valeurs NFIRE du registre (W0 = 3,25 cm, C_r = π)
c3 = copy.deepcopy(cfg); c3.w0_gs_m = 0.0325; c3.c_r = np.pi
for z, sb, si in ((0, 0.035, 0.17), (60, 0.165, 0.46)):
    st = n1.uplink_stats(np.radians(z), 550e3, c3)
    check(f"A3 registre P1-C3-CERT ζ={z}°", abs(st["sigma_b2"] / sb - 1) < 0.03 and abs(st["sigma_i2"] / si - 1) < 0.03,
          f"σ²_B={st['sigma_b2']:.3f} σ²_I={st['sigma_i2']:.3f}")
# A4 : errance angulaire Fig. 2 (≈ 7 et 10 µrad pour W0 ≈ 3 cm)
for z, ang in ((0, 7.1), (60, 10.0)):
    st = n1.uplink_stats(np.radians(z), 550e3, c3)
    check(f"A4 errance angulaire ζ={z}°", abs(np.sqrt(st["rc2"]) / st["L"] * 1e6 - ang) < 0.3)
# A5 : C6 = limite Λ = 0, Θ̄ = 1 de l'Éq. 33
k = cfg.k; H = 550e3
f_sph = lambda h: n1.v19.cn2_hv(h, cfg.v19) * np.real((1 - h / H) ** (5 / 6) * (1j * (h / H)) ** (5 / 6))
mu = n1.quad(f_sph, 0, 30e3, limit=400, points=[100, 1000, 10000])[0]
c6v = 2.25 * k ** (7 / 6) * n1.quad(lambda h: n1.v19.cn2_hv(h, cfg.v19) * (h * (H - h) / H) ** (5 / 6),
                                     0, 30e3, limit=400, points=[100, 1000, 10000])[0]
check("A5 Éq. 33 avec Λ=0, Θ̄=1 ⇒ C6 (onde sphérique)", abs(8.70 * mu * k ** (7 / 6) * H ** (5 / 6) / c6v - 1) < 1e-3)

print("B — Récepteur et MCS")
loss = n1.implementation_loss_db(cfg)
n_ideal = n1.erfcinv(2e-8) ** 2 / 2
check("B1 perte = 10·log10(48 / n_idéal), n_idéal = 7,874 ph/bit ⇒ 7,85 dB (registre : 7,7 dB arrondi, à corriger)",
      abs(n_ideal - 7.874) < 0.001 and abs(loss - 10 * np.log10(48 / n_ideal)) < 1e-12, f"{loss:.3f} dB")
p_meas = cfg.ppb_measured * cfg.photon_j * 1.8e9
check("B2 48 photons/bit à 1,8 Gbit/s ⇔ −47,9 dBm (ICSOS 2012 : −47,8)",
      abs(10 * np.log10(p_meas * 1e3) + 47.85) < 0.2)
esn0 = n1.esn0_from_power(p_meas, 1.8e9, cfg)
ber = 0.5 * n1.erfc(np.sqrt(esn0))
check("B3 E_s/N0 au point mesuré ⇒ BER BPSK = 1e-8", abs(ber / 1e-8 - 1) < 1e-6, f"{ber:.3e}")
thr = n1.bpsk_thresholds_lin()
check("B4 seuils strictement croissants, 11 rendements", len(thr) == 11 and np.all(np.diff(thr) > 0))
check("B5 BPSK = QPSK − 3,01 dB", np.allclose(10 * np.log10(thr) - np.array(n1.QPSK_ESN0_DB), -3.0103, atol=1e-3))

print("C — Échantillonnage")
rng = np.random.default_rng(3)
x = n1.sample_lognormal(0.3, 400_000, rng)
check("C1 log-normale : moyenne 1, variance σ²_I", abs(x.mean() - 1) < 0.005 and abs(x.var() / 0.3 - 1) < 0.05)
check("C2 captation gaussienne ≤ puissance émise",
      0 < n1.gaussian_rx_power(1.0, 5.0, 0.124) < 1.0)

print("D — Politiques V1.3.3 (invariants, 5 époques, aucune valeur de J publiée)")
small = copy.deepcopy(cfg); small.episode_epochs = 5; small.n_mc = 500
ch = n1.N1Channel(small, seed=1)
eps = [ch.epoch(t * small.dt_s) for t in range(5)]
dm = n1.DecisionModel(n1.ChannelModel(small))
check("D1 visibilité limitée à ζ ≤ 60°", all(np.all(np.degrees(e["zen"][e["vis"]]) <= 60.0) for e in eps))
ok = True
for ep in eps:
    allg = [n1.realized_goodput(ep, i, m, small) for i in range(small.n_relays) for m in range(-1, 11)]
    res1 = n1.evaluate_trajectory(ch, [ep], dm)
    ok &= res1["clairvoyante"]["J"] >= max(allg) - 1e-12
check("D2 clairvoyante = meilleure action réalisée à chaque époque", ok)
# D3 : absence de consultation du canal courant — Q identique pour deux réalisations de même géométrie
ch_a, ch_b = n1.N1Channel(small, seed=1), n1.N1Channel(small, seed=1)
ea = ch_a.epoch(0.0)
ch_b.rng = np.random.default_rng(999)          # même géométrie, évanouissements différents
eb = ch_b.epoch(0.0)
check("D3 Q_t ne dépend que de la géométrie (pas du canal réalisé)",
      not np.array_equal(ea["g1"], eb["g1"]) and np.array_equal(dm.q_relays(ea), dm.q_relays(eb)))
P = dm.p_hop1(np.radians(30.0))
check("D4 P(γ̄₁ ≥ seuil) décroissante avec le rendement", np.all(np.diff(P) <= 1e-12))
res = n1.evaluate_trajectory(ch, eps, dm)
check("D5 six comparateurs sur les mêmes réalisations",
      set(res) == {"causale", "clairvoyante", "r_fixe", "m_fixe", "statique", "directe"})
check("D6 clairvoyante ≥ chaque comparateur relayé en J",
      all(res["clairvoyante"]["J"] >= res[p]["J"] - 1e-12 for p in res if p != "directe"))
dm2 = n1.DecisionModel(n1.ChannelModel(small))
check("D7 modèle de décision reproductible (flux MC à graine fixe)",
      np.array_equal(dm.p_hop1(np.radians(30.0)), dm2.p_hop1(np.radians(30.0))))

print("E — Critères et gouvernance")
rng = np.random.default_rng(0)
st = n1.contrast_stats(1.1 + 0.01 * rng.normal(size=20), np.ones(20))
check("E1 contraste : écart relatif, IC BCa, p inversion de signe", abs(st["rel"] - 0.1) < 0.01 and st["ci"][0] > 0 and st["p"] < 0.01)
check("E2 Holm step-down", list(n1.holm([0.01, 0.04, 0.03, 0.2])) == [True, False, False, False]
      and list(n1.holm([0.001, 0.012, 0.02, 0.04])) == [True, True, True, True])
try:
    n1.run_grid(); ok = False
except PermissionError:
    ok = True
check("E3 run_grid refuse sans jeton d'autorisation", ok)
check("E4 grille = 4 configurations (W0 × C_r), 20 graines", len(n1.GRID) == 4 and len(n1.SEEDS) == 20)
check("E5 V19 importée sous empreinte gelée", n1.V19_SHA256.startswith("a8e6dcbb"))
c77 = copy.deepcopy(cfg); c77.impl_loss_db = 7.7
check("E6 perte d'implémentation paramétrable (point ouvert A-1 : 7,7 ou 7,85 dB)",
      n1.implementation_loss_db(c77) == 7.7 and abs(n1.implementation_loss_db(cfg) - 7.85) < 0.01)

T, F = True, False
check("E7 classification : A (ΔJ_R et ΔJ_tot sur ≥ 2 configurations communes)",
      n1.classify([T, T, F, F], [F]*4, [T, T, F, F], [T]*4) == "A")
check("E8 classification : ΔJ_R sur 2 configurations, ΔJ_tot nulle part ⇒ C (étape 5)",
      n1.classify([T, T, F, F], [F]*4, [F]*4, [T]*4) == "C")
check("E8b classification : ΔJ_R et ΔJ_tot sur des configurations disjointes ⇒ Indéterminée",
      n1.classify([T, T, F, F], [F]*4, [F, F, T, F], [T]*4) == "Indéterminée")
check("E9 classification : A partielle, A non exploitable, B, C",
      n1.classify([T, F, F, F], [F]*4, [T]*4, [T]*4) == "A partielle"
      and n1.classify([F]*4, [F]*4, [T]*4, [T, F, F, F]) == "A non exploitable"
      and n1.classify([F]*4, [F, T, F, F], [T]*4, [F]*4) == "B"
      and n1.classify([F]*4, [F]*4, [F]*4, [F]*4) == "C")
fam = [n1.contrast_stats(1.1 + 0.01 * rng.normal(size=20), np.ones(20)) for _ in range(3)] + \
      [n1.contrast_stats(1.02 + 0.01 * rng.normal(size=20), np.ones(20))]
check("E10 significativité = Holm ∧ ≥ 5 % ∧ BCa > 0 (un gain de 2 % est rejeté)",
      list(n1.significance_family(fam)) == [True, True, True, False])

# E11–E13 : contrôles analytiques du test exact (aucune grille).
check("E11 20 différences strictement positives : p exacte = 2^-20",
      n1.sign_flip_pvalue_exact(np.ones(20)) == 1 / (2 ** 20))
check("E12 deux différences opposées : p unilatérale = 3/4",
      n1.sign_flip_pvalue_exact(np.array([1.0, -1.0])) == 0.75)
try:
    n1.sign_flip_pvalue_exact(np.ones(21))
    rejected = False
except ValueError:
    rejected = True
check("E13 refus des échantillons dépassant 20 paires", rejected)

# E14–E16 : clause d'asymétrie (V1.3.3)
c_sym = n1.contrast_stats(1.1 + 0.01 * np.random.default_rng(5).normal(size=20), np.ones(20))
d_sk = np.r_[np.full(19, 0.12), 3.0]
c_sk = n1.contrast_stats(1.0 + d_sk, np.ones(20))
check("E14 asymétrie calculée et retournée par contrast_stats", abs(c_sym["skew"]) < 1 and c_sk["skew"] > 1,
      f"{c_sym['skew']:.2f} / {c_sk['skew']:.2f}")
fl = n1.asymmetry_flags([c_sym, c_sk, c_sk], np.array([True, True, False]))
check("E15 marquage « à confirmer » seulement si significatif et |asymétrie| > 1", list(fl) == [False, True, False])
# E16 : boucle de campagne en miniature (1 config, 2 graines, 3 époques) — structure seulement
import tempfile, json, os
tmp = tempfile.mkdtemp()
mini = copy.deepcopy(cfg); mini.n_mc = 200
summ = n1._campaign(n1.GRID[:1], [0, 1], 3, tmp, base_cfg=mini)
files_ok = all(os.path.exists(os.path.join(tmp, f)) for f in ("runs.csv", "contrasts.csv", "summary.json", "SHA256SUMS.txt"))
check("E16 boucle de campagne : 4 fichiers produits, issue et empreintes présentes (miniature, valeurs non lues)",
      files_ok and "issue" in summ and summ["code_sha256"] and summ["v19_sha256"] == n1.V19_SHA256)
try:
    n1.run_grid("jeton-invalide"); ok = False
except PermissionError:
    ok = True
check("E17 run_grid refuse un jeton invalide", ok)

print(f"\n{len(PASSED)} vérifications réussies, 0 échec.")
