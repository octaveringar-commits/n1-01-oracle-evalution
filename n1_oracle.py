"""
n1_oracle.py — Oracle exploratoire ΔJ, fiche N1-01 V1.3.3 (registre V18.4.3).
=============================================================================

STATUT : préparé pour AUDIT STATIQUE. Exécution de la grille NON AUTORISÉE
(N1-01 V1.3.3, §10). run_grid() refuse de s'exécuter sans le
jeton d'autorisation explicite ; seuls les tests unitaires peuvent tourner.

Niveau G-02 : N1 (exploratoire, non publiable). V19 révision 2 est importée
SANS modification (géométrie, profil HV5/7) ; empreinte vérifiée à l'import.

Architecture (décision P1-C-DF-01) : Sol → LEO (DF) → GEO, 1064 nm,
BPSK homodyne cohérent, action a_t = (R_t, M_t), MCS commun aux deux sauts.

Traçabilité : [S] source, [D] déduction, [H] hypothèse, [ND] non documenté.

Sources
-------
[A06] Andrews, Phillips, Sasiela, Parenti, Opt. Eng. 45(7) 076001, 2006,
      DOI 10.1117/1.2219470 (équations montantes certifiées, P1-C3-CERT).
[AL22] Alaluf, Perdigues Armengol, CEAS Space J. 14:227-238, 2022,
      DOI 10.1007/s12567-021-00392-2 (cas non suivi, jitter OGS).
[F11] Fields et al., Proc. SPIE 8184, 81840D, 2011 (NFIRE : 6,5 / 12,4 cm).
[S08] Smutny et al., Proc. SPIE 6877, 2008 (LCT : 0,7 / 2,2 W, 135 mm).
[I12] ICSOS 2012, LCT Alphasat : 48 photons/bit à BER 1e-8, 1,8 Gbit/s.
[ETSI] ETSI EN 302 307-1 V1.4.1 (2014-11), clause 6, Table 13 (p. 36) — seuils QPSK vérifiés
      sur le PDF (SHA-256 19077f80…d8475d), PER QEF = 1e-7, AWGN, synchronisation parfaite.
"""
from __future__ import annotations

import copy
import hashlib
import os
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import quad
from scipy.special import erfcinv, erfc

_HERE = os.path.dirname(os.path.abspath(__file__))
_V19 = os.path.join(os.path.dirname(_HERE), "v19")
sys.path.insert(0, _V19)
V19_SHA256 = "a8e6dcbb77b0a21b40530cf285b455b667dc8c3e9f6ba685de69376135c3f851"
with open(os.path.join(_V19, "fso_env.py"), "rb") as _f:
    if hashlib.sha256(_f.read()).hexdigest() != V19_SHA256:
        raise RuntimeError("fso_env.py (V19 rév. 2) modifié : empreinte différente du gel")
import fso_env as v19  # noqa: E402

H_PLANCK = 6.62607015e-34
C_LIGHT = 299_792_458.0
AUTH_TOKEN = "N1-01-V1.3.3-EXECUTION-AUTORISEE"   # à fournir explicitement par l'auteur

# Seuils E_s/N0 idéaux (dB) QPSK DVB-S2, FECFRAME 64 800, PER QEF 1e-7, AWGN [ETSI, Table 13]
# A-4 clos : 11 valeurs vérifiées le 20/10/2026 sur ETSI EN 302 307-1 V1.4.1, p. 36.
ETSI_PDF_SHA256 = "19077f80420eb17519fcf6f83e2847d7a877940b134c5306814b8c6f60d8475d"
LDPC_RATES = (1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, 9/10)
QPSK_ESN0_DB = (-2.35, -1.24, -0.30, 1.00, 2.23, 3.10, 4.03, 4.68, 5.18, 6.20, 6.42)
BPSK_OFFSET_DB = -10 * np.log10(2.0)   # [D] BPSK = une voie de QPSK : E_s/N0 requis −3,01 dB


@dataclass
class N1Config:
    # ── Saut 1 : sol → LEO (NFIRE) ─────────────────────────────────
    lambda_m: float = 1064e-9           # [S] F11, S08
    pt_gs_w: float = 0.7                # [S] S08 puissance LCT ; [H] même valeur au sol
    d_tx_gs_m: float = 0.065            # [S] F11 ouverture sol ; [H] émission = réception
    w0_gs_m: float = 0.065 / 2 ** 1.5   # décision 4 : W0 ∈ {D/2^{3/2}, D/2} [H]
    d_rx_leo_m: float = 0.124           # [S] F11 (ambiguïté 124/125/135 mm déclarée)
    rs1_baud: float = 5.625e9           # [H] BPSK aux débits mesurés (décision N1-01)
    jitter_gs_rad: float = 2e-6         # décision 4 ; [S] AL22 : 0,7–1,7 µrad (OGS 1 m)
    c_r: float = np.pi                  # décision 4 : C_r ∈ {π ; 3,86} [S] A06
    tau_zenith: float = 0.10            # [H] H-ATM, hérité V19
    # ── Saut 2 : LEO → GEO (EDRS) ──────────────────────────────────
    pt_leo_w: float = 2.2               # [S] S08
    d_tx_leo_m: float = 0.135           # [S] S08
    d_rx_geo_m: float = 0.135           # [S] S08
    rs2_baud: float = 1.8e9             # [H] débit utile EDRS en BPSK non codé
    jitter_isl_rad: float = 1e-6        # [H] H-PE2, hérité V19
    # ── Récepteur ──────────────────────────────────────────────────
    ppb_measured: float = 48.0          # [S] I12 photons/bit à BER 1e-8
    ber_measured: float = 1e-8          # [S] I12
    impl_loss_db: float | None = None   # A-1 tranché : None ⇒ 7,85 dB dérivé de I12, η_eff = 1 ;
                                        # 7,7 dB conservé comme option de sensibilité seulement
    # ── Géométrie, statistique ─────────────────────────────────────
    zenith_max_deg: float = 60.0        # décision 4 : validité A06 (ζ ≤ 60°)
    n_relays: int = 10                  # T4
    dt_s: float = 0.1                   # T4
    n_subslots: int = 64                # H-BLK, hérité V19
    episode_epochs: int = 1000          # V1.3.3 §9
    n_mc: int = 4000                    # tirages Monte-Carlo du modèle pour Q_t (§7), flux séparé
    mc_seed: int = 20261020             # graine du modèle de décision, indépendante des réalisations
    h_turb_top_m: float = 30e3
    v19: v19.Config = field(default_factory=lambda: copy.deepcopy(v19.CONFIG))

    @property
    def k(self) -> float:
        return 2 * np.pi / self.lambda_m

    @property
    def photon_j(self) -> float:
        return H_PLANCK * C_LIGHT / self.lambda_m


# ─────────────────────────────────────────────────────────────────────
# Récepteur et MCS
# ─────────────────────────────────────────────────────────────────────
def implementation_loss_db(cfg: N1Config) -> float:
    """[D] Perte vs BPSK homodyne idéal : BER = ½ erfc(√(2 n)), n photons/bit."""
    if cfg.impl_loss_db is not None:
        return float(cfg.impl_loss_db)
    n_ideal = erfcinv(2 * cfg.ber_measured) ** 2 / 2.0
    return 10 * np.log10(cfg.ppb_measured / n_ideal)


def esn0_from_power(p_w, rs_baud: float, cfg: N1Config):
    """E_s/N0 (linéaire) = 2·n_s / L_impl, n_s = P/(hν R_s) photons par symbole [D]."""
    n_s = np.asarray(p_w) / (cfg.photon_j * rs_baud)
    return 2.0 * n_s / 10 ** (implementation_loss_db(cfg) / 10)


def bpsk_thresholds_lin() -> np.ndarray:
    return 10 ** ((np.asarray(QPSK_ESN0_DB) + BPSK_OFFSET_DB) / 10)


# ─────────────────────────────────────────────────────────────────────
# Canal montant certifié [A06]
# ─────────────────────────────────────────────────────────────────────
def mu0(cfg: N1Config) -> float:
    return quad(lambda h: v19.cn2_hv(h, cfg.v19), 0, cfg.h_turb_top_m,
                limit=400, points=[100, 1000, 10000])[0]


def uplink_stats(zenith_rad: float, H: float, cfg: N1Config) -> dict:
    """
    Faisceau collimaté non suivi [A06], h0 = 0, géométrie plane (sec ζ) :
      r0 (Éq. 22), ⟨r_c²⟩ (Éq. 20), σ²_pe (Éq. 30), W (Éq. 3),
      μ_u1 (Éq. 33), σ²_B (Éq. 37), σ²_I,l non suivi forme forte (Éq. 44).
    """
    k, W0 = cfg.k, cfg.w0_gs_m
    sec = 1.0 / np.cos(zenith_rad)
    L = H * sec
    r0 = (0.423 * mu0(cfg) * k ** 2 * sec) ** (-3 / 5)
    rc2 = 0.54 * H ** 2 * sec ** 2 * (cfg.lambda_m / (2 * W0)) ** 2 * (2 * W0 / r0) ** (5 / 3)
    x = cfg.c_r ** 2 * W0 ** 2 / r0 ** 2
    spe2 = rc2 * (1 - (x / (1 + x)) ** (1 / 6))
    lam0 = 2 * L / (k * W0 ** 2)
    theta, lam = 1 / (1 + lam0 ** 2), lam0 / (1 + lam0 ** 2)
    thb = 1 - theta
    W = W0 * np.sqrt(1 + lam0 ** 2)

    def f(h):
        xi = 1 - h / H
        return v19.cn2_hv(h, cfg.v19) * np.real(
            xi ** (5 / 6) * (lam * xi + 1j * (1 - thb * xi)) ** (5 / 6) - lam ** (5 / 6) * xi ** (5 / 3))

    mu_u1 = quad(f, 0, cfg.h_turb_top_m, limit=400, points=[100, 1000, 10000])[0]
    sb2 = 8.70 * mu_u1 * k ** (7 / 6) * H ** (5 / 6) * sec ** (11 / 6)
    a_pe = np.sqrt(spe2) / L
    wander = 5.95 * H ** 2 * sec ** 2 * (2 * W0 / r0) ** (5 / 3) * (a_pe / W) ** 2
    strong = np.exp(0.49 * sb2 / (1 + 0.56 * (1 + theta) * sb2 ** 1.2) ** (7 / 6)
                    + 0.51 * sb2 / (1 + 0.69 * sb2 ** 1.2) ** (5 / 6)) - 1
    # W_LT² = W_ST² + ⟨r_c²⟩ (Éq. 13) avec W_ST ≈ W : étalement petite échelle négligé
    # [H-WST], justifié par W0/r0 ≤ 0,4 (régime dominé par l'inclinaison, A06 Fig. 3).
    return dict(L=L, r0=r0, rc2=rc2, spe2=spe2, W=W, Lambda=lam, sigma_b2=sb2,
                sigma_i2=wander + strong, w_lt=np.sqrt(W ** 2 + rc2))


def gaussian_rx_power(pt_w, w_m, d_rx_m):
    """Puissance captée par une ouverture centrée, faisceau gaussien de rayon w (1/e²)."""
    a = d_rx_m / 2
    return pt_w * (1 - np.exp(-2 * a ** 2 / w_m ** 2))


def sample_lognormal(sigma_i2, size, rng):
    """[S] A06 Fig. 3 : statistique log-normale de moyenne 1 pour W0/r0 ≪ 1."""
    s2 = np.log1p(sigma_i2)
    return np.exp(rng.normal(-s2 / 2, np.sqrt(s2), size))


# ─────────────────────────────────────────────────────────────────────
# Modèle de canal (partagé par le simulateur ET par le modèle de décision)
# ─────────────────────────────────────────────────────────────────────
class ChannelModel:
    """Fonctions d'échantillonnage du canal ; un seul code pour les réalisations et pour Q_t."""

    def __init__(self, cfg: N1Config):
        self.cfg = cfg
        self._up_cache: dict = {}

    def up(self, zenith_rad: float, H: float):
        # [H-ARR] Arrondi de cache à 0,5° ; erreur décisionnelle non quantifiée.
        key = (round(np.degrees(zenith_rad) * 2) / 2, H)
        if key not in self._up_cache:
            self._up_cache[key] = uplink_stats(np.radians(key[0]), H, self.cfg)
        return self._up_cache[key]

    def uplink_esn0(self, zenith_rad, H, rs_baud, d_rx_m, size, rng):
        """Montée : P_ref (W_LT, atmosphère) × jitter mécanique × log-normale (σ²_I non suivi)."""
        cfg = self.cfg
        st = self.up(zenith_rad, H)
        latm = np.exp(-cfg.tau_zenith / np.cos(zenith_rad))
        p = gaussian_rx_power(cfg.pt_gs_w * latm, st["w_lt"], d_rx_m)
        rad = rng.rayleigh(cfg.jitter_gs_rad * st["L"], size)
        p = p * np.exp(-2 * rad ** 2 / st["w_lt"] ** 2) * sample_lognormal(st["sigma_i2"], size, rng)
        return esn0_from_power(p, rs_baud, cfg)

    def isl_esn0(self, d2_m, size, rng):
        """LEO → GEO dans le vide : diffraction (W0 = D/2, convention V19) × jitter ISL."""
        cfg = self.cfg
        w2 = cfg.lambda_m * d2_m / (np.pi * cfg.d_tx_leo_m / 2)
        rad = rng.rayleigh(cfg.jitter_isl_rad * d2_m, size)
        p = gaussian_rx_power(cfg.pt_leo_w, w2, cfg.d_rx_geo_m) * np.exp(-2 * rad ** 2 / w2 ** 2)
        return esn0_from_power(p, cfg.rs2_baud, cfg)


class DecisionModel:
    """
    Q_t (V1.3.3 §7) : probabilité que la moyenne d'époque franchisse chaque seuil,
    estimée par Monte-Carlo du MODÈLE, sur un flux aléatoire distinct des réalisations.
    Ne dépend que de la géométrie courante (ζ, d₂) et des paramètres Θ : aucun accès
    au canal courant ni futur.
    """

    def __init__(self, cm: ChannelModel):
        self.cm, self.cfg = cm, cm.cfg
        self.thr = bpsk_thresholds_lin()
        self._p1: dict = {}
        self._p2: dict = {}
        self._pd = None

    def _tail(self, means):
        means = np.sort(means)
        return 1.0 - np.searchsorted(means, self.thr, side="left") / len(means)

    def p_hop1(self, zenith_rad):
        # [H-ARR] Approximation angulaire à 0,5° ; effet sur Q non quantifié.
        key = round(np.degrees(zenith_rad) * 2) / 2
        if key not in self._p1:
            rng = np.random.default_rng([self.cfg.mc_seed, 1, int(key * 2)])
            g = self.cm.uplink_esn0(np.radians(key), self.cfg.v19.h_leo_m, self.cfg.rs1_baud,
                                    self.cfg.d_rx_leo_m, (self.cfg.n_mc, self.cfg.n_subslots), rng)
            self._p1[key] = self._tail(g.mean(axis=1))
        return self._p1[key]

    def p_hop2(self, d2_m):
        key = int(round(d2_m / 1e5))  # [H-ARR] pas de 100 km ; effet sur Q non quantifié
        if key not in self._p2:
            rng = np.random.default_rng([self.cfg.mc_seed, 2, key])
            g = self.cm.isl_esn0(key * 1e5, (self.cfg.n_mc, self.cfg.n_subslots), rng)
            self._p2[key] = self._tail(g.mean(axis=1))
        return self._p2[key]

    def p_direct(self, zenith_geo):
        if self._pd is None:
            rng = np.random.default_rng([self.cfg.mc_seed, 3])
            g = self.cm.uplink_esn0(zenith_geo, self.cfg.v19.h_geo_m, self.cfg.rs2_baud,
                                    self.cfg.d_rx_geo_m, (self.cfg.n_mc, self.cfg.n_subslots), rng)
            self._pd = self._tail(g.mean(axis=1))
        return self._pd

    def q_relays(self, geom):
        """Q[i, m] = 1,8·r_m·P₁·P₂ pour les relais admissibles, 0 sinon (§7)."""
        rates = np.asarray(LDPC_RATES) * min(self.cfg.rs1_baud, self.cfg.rs2_baud) / 1e9
        Q = np.zeros((self.cfg.n_relays, len(LDPC_RATES)))
        for i in np.where(geom["vis"])[0]:
            Q[i] = rates * self.p_hop1(geom["zen"][i]) * self.p_hop2(geom["d2"][i])
        return Q

    def q_direct(self, geom):
        return np.asarray(LDPC_RATES) * self.cfg.rs2_baud / 1e9 * self.p_direct(geom["zen_geo"])


# ─────────────────────────────────────────────────────────────────────
# Simulateur de trajectoire (réalisations)
# ─────────────────────────────────────────────────────────────────────
class N1Channel:
    def __init__(self, cfg: N1Config, seed: int):
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        vc = copy.deepcopy(cfg.v19)
        vc.n_relays = cfg.n_relays
        self.vc = vc
        self.gs, self.geo = v19.gs_position(vc), v19.geo_position(vc)
        self.const = v19.Constellation(vc, self.rng)
        self.cm = ChannelModel(cfg)

    def geometry(self, t_s: float):
        """𝒢_t : géométrie seule (connue de la politique, éphémérides)."""
        sats = self.const.positions(t_s)
        zen = v19.zenith_angle(self.gs, sats)
        vis = (np.degrees(zen) <= self.cfg.zenith_max_deg) & v19.line_of_sight_clear(self.gs, sats) \
            & v19.line_of_sight_clear(sats, self.geo)
        d2 = np.linalg.norm(self.geo - sats, axis=-1)
        zen_geo = v19.zenith_angle(self.gs, self.geo[None])[0]
        return dict(zen=zen, vis=vis, d2=d2, zen_geo=zen_geo, sats=sats)

    def epoch(self, t_s: float):
        """Réalisation : E_s/N0 par sous-créneau sur les deux sauts et la liaison directe."""
        cfg, K, rng = self.cfg, self.cfg.n_subslots, self.rng
        g = self.geometry(t_s)
        N = cfg.n_relays
        g1, g2 = np.zeros((N, K)), np.zeros((N, K))
        for i in np.where(g["vis"])[0]:
            g1[i] = self.cm.uplink_esn0(g["zen"][i], cfg.v19.h_leo_m, cfg.rs1_baud, cfg.d_rx_leo_m, K, rng)
            g2[i] = self.cm.isl_esn0(g["d2"][i], K, rng)
        gd = self.cm.uplink_esn0(g["zen_geo"], cfg.v19.h_geo_m, cfg.rs2_baud, cfg.d_rx_geo_m, K, rng)
        g.update(g1=g1, g2=g2, gd=gd)
        return g


# ─────────────────────────────────────────────────────────────────────
# Décodage et politiques (V1.3.3 §6–7)
# ─────────────────────────────────────────────────────────────────────
def best_rate_index(esn0_mean_lin: float, thr: np.ndarray) -> int:
    ok = np.where(esn0_mean_lin >= thr)[0]
    return int(ok[-1]) if len(ok) else -1


def realized_goodput(ep, relay, m, cfg, direct=False) -> float:
    """G_t (§6) : décodage sur la moyenne linéaire de l'époque [H], deux sauts (DF, MCS commun)."""
    thr = bpsk_thresholds_lin()
    if m is None or m < 0:
        return 0.0
    if direct:
        return LDPC_RATES[m] * cfg.rs2_baud / 1e9 if ep["gd"].mean() >= thr[m] else 0.0
    if relay is None or not ep["vis"][relay]:
        return 0.0
    ok = ep["g1"][relay].mean() >= thr[m] and ep["g2"][relay].mean() >= thr[m]
    return LDPC_RATES[m] * min(cfg.rs1_baud, cfg.rs2_baud) / 1e9 if ok else 0.0


def highest_elevation_relay(geom):
    if not geom["vis"].any():
        return None
    return int(np.argmin(np.where(geom["vis"], geom["zen"], np.inf)))


def argmax_q(Q):
    """(R, m) maximisant Q ; (None, −1) si Q ≡ 0."""
    if Q.max() <= 0:
        return None, -1
    i, m = np.unravel_index(int(np.argmax(Q)), Q.shape)
    return int(i), int(m)


def evaluate_trajectory(ch: N1Channel, epochs, dm: DecisionModel) -> dict:
    """
    Six comparateurs sur les MÊMES réalisations (§7). Les politiques causales ne lisent que
    la géométrie 𝒢_t et le modèle Θ (via dm) ; seule la clairvoyante lit la réalisation.
    Rendement constant (M-fixe, statique) : choisi EX ANTE par maximisation de Σ_t Q_t
    (éphémérides + modèle), sans consulter les réalisations (§9, absence de consultation du futur).
    """
    cfg = ch.cfg
    T = len(epochs)
    Qs = [dm.q_relays(ep) for ep in epochs]
    Qd = [dm.q_direct(ep) for ep in epochs]
    he = [highest_elevation_relay(ep) for ep in epochs]
    m_mfix = int(np.argmax(sum(Q.max(axis=0) for Q in Qs)))
    m_stat = int(np.argmax(sum((Q[r] if r is not None else np.zeros(Q.shape[1])) for Q, r in zip(Qs, he))))
    names = ("causale", "clairvoyante", "r_fixe", "m_fixe", "statique", "directe")
    out = {p: np.zeros(T) for p in names}
    acts_hist = {p: [] for p in names}
    for t, ep in enumerate(epochs):
        Q = Qs[t]
        r_he = he[t]
        # clairvoyante : meilleure action RÉALISÉE (borne supérieure)
        best = (None, -1, 0.0)
        for i in np.where(ep["vis"])[0]:
            for m in range(len(LDPC_RATES) - 1, -1, -1):
                g = realized_goodput(ep, int(i), m, cfg)
                if g > best[2]:
                    best = (int(i), m, g)
                    break
        acts = {
            "causale": argmax_q(Q),
            "clairvoyante": best[:2],
            "r_fixe": (r_he, int(np.argmax(Q[r_he])) if r_he is not None and Q[r_he].max() > 0 else -1),
            "m_fixe": ((int(np.argmax(Q[:, m_mfix])) if Q[:, m_mfix].max() > 0 else None), m_mfix),
            "statique": (r_he, m_stat),
            "directe": (None, int(np.argmax(Qd[t])) if Qd[t].max() > 0 else -1),
        }
        for p, (r, m) in acts.items():
            out[p][t] = realized_goodput(ep, r, m, cfg, direct=(p == "directe"))
            acts_hist[p].append((r, m))
    res = {}
    for p in names:
        h = acts_hist[p]
        sr = sum(a[0] != b[0] for a, b in zip(h[1:], h[:-1])) / max(T - 1, 1)
        sm = sum(a[1] != b[1] for a, b in zip(h[1:], h[:-1])) / max(T - 1, 1)
        res[p] = dict(J=float(out[p].mean()), outage=float((out[p] == 0).mean()),
                      switch_relay=sr, switch_mcs=sm)
    return res


def simulate(cfg: N1Config, seed: int) -> dict:
    ch = N1Channel(cfg, seed)
    dm = DecisionModel(ChannelModel(cfg))
    epochs = [ch.epoch(t * cfg.dt_s) for t in range(cfg.episode_epochs)]
    return evaluate_trajectory(ch, epochs, dm)


# ─────────────────────────────────────────────────────────────────────
# Grille préenregistrée et critères (V1.3.3 §3, §8)
# ─────────────────────────────────────────────────────────────────────
GRID = [dict(w0_gs_m=w0, c_r=cr) for w0 in (0.065 / 2 ** 1.5, 0.065 / 2) for cr in (np.pi, 3.86)]
SEEDS = list(range(20))
THRESHOLD_REL = 0.05
ALPHA = 0.05


def sign_flip_pvalue_exact(d, chunk_size=65536) -> float:
    """[H-STAT] Test unilatéral par inversion exhaustive des signes appariés.

    Statistique : moyenne des différences d = a - b ; alternative positive.
    p = #{s ∈ {-1,+1}^n : mean(s*d) >= mean(d)} / 2^n.
    Exact sous invariance de la loi jointe de d par inversions de signes
    (notamment différences indépendantes et symétriques sous H0).
    Cette condition est PLUS FORTE que la seule hypothèse E[d] <= 0.
    Le test n'est pas certifié exact pour H0 portant seulement sur la moyenne.
    La campagne préenregistrée impose n=20 ; la limite protège l'énumération.
    """
    d = np.asarray(d, dtype=float)
    if d.ndim != 1 or not 1 <= len(d) <= 20 or not np.all(np.isfinite(d)):
        raise ValueError("Différences finies appariées, 1 <= n <= 20 requises")
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError("chunk_size doit être un entier strictement positif")
    n = len(d)
    total = 1 << n
    observed = float(d.sum())
    extreme = 0
    bit_positions = np.arange(n, dtype=np.uint64)
    for start in range(0, total, chunk_size):
        masks = np.arange(start, min(start + chunk_size, total), dtype=np.uint64)
        flipped = ((masks[:, None] >> bit_positions) & np.uint64(1)).astype(bool)
        sums = np.where(flipped, -d, d).sum(axis=1)
        extreme += int(np.count_nonzero(sums >= observed - 1e-12 * max(1.0, abs(observed))))
    return extreme / total


def contrast_stats(a, b, n_boot=10000, seed=42):
    """Écart relatif, IC BCa 95 % et p unilatérale exacte par inversion de signe.

    [H-STAT] Exactitude conditionnelle à l'invariance des signes sous H0.
    Le BCa décrit l'incertitude sur la moyenne ; il ne valide pas l'hypothèse
    d'échangeabilité requise par le test de permutation.
    """
    from scipy.stats import bootstrap
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("a et b doivent être deux vecteurs appariés de même longueur")
    d = a - b
    p = sign_flip_pvalue_exact(d)
    ref = float(np.mean(b))
    rel = float(np.mean(d)) / ref if ref > 0 else np.inf
    ci = bootstrap((d,), np.mean, method="BCa", confidence_level=0.95,
                   n_resamples=n_boot, random_state=seed).confidence_interval
    # V1.3.3 — clause d'asymétrie : coefficient d'asymétrie corrigé (Fisher-Pearson ajusté)
    from scipy.stats import skew
    sk = float(skew(d, bias=False)) if np.ptp(d) > 0 else 0.0
    return dict(rel=rel, ci=(float(ci.low), float(ci.high)), p=p, skew=sk)


def holm(pvals, alpha=ALPHA):
    """Rejets de Holm–Bonferroni (step-down) dans l'ordre d'entrée."""
    p = np.asarray(pvals)
    order = np.argsort(p)
    rej = np.zeros(len(p), bool)
    for k, idx in enumerate(order):
        if p[idx] <= alpha / (len(p) - k):
            rej[idx] = True
        else:
            break
    return rej


def significance_family(per_config: list[dict]) -> np.ndarray:
    """
    A-3 (validé) : famille = les 4 configurations d'UN contraste.
    Significatif ⇔ rejet de Holm à 5 % ET écart relatif ≥ 5 % ET borne basse BCa > 0.
    per_config : liste de sorties de contrast_stats, dans l'ordre de GRID.
    """
    rej = holm([c["p"] for c in per_config])
    return np.array([r and c["rel"] >= THRESHOLD_REL and c["ci"][0] > 0
                     for r, c in zip(rej, per_config)])


SKEW_LIMIT = 1.0


def asymmetry_flags(per_config: list[dict], sig: np.ndarray) -> np.ndarray:
    """
    V1.3.3 — clause d'asymétrie (préenregistrée) : un écart SIGNIFICATIF est marqué
    « à confirmer » si |asymétrie des 20 écarts appariés| > 1. Le marquage ne retire
    pas la significativité ; il est reporté avec l'issue.
    """
    return np.array([bool(s) and abs(c["skew"]) > SKEW_LIMIT for s, c in zip(sig, per_config)])


def classify(sig_R, sig_M, sig_tot, sig_clair_vs_R) -> str:
    """
    A-2 : règle ordonnée sur 4 configurations (première condition vraie retenue).
    Les arguments sont des vecteurs booléens de significativité (sortie de significance_family).
    """
    sig_R, sig_M, sig_tot, sig_c = map(np.asarray, (sig_R, sig_M, sig_tot, sig_clair_vs_R))
    if np.sum(sig_R & sig_tot) >= 2:
        return "A"
    if np.sum(sig_R) == 1:
        return "A partielle"
    if not sig_R.any() and sig_c.any():
        return "A non exploitable"
    if not sig_R.any() and sig_M.any():
        return "B"
    if not sig_tot.any():
        return "C"
    return "Indéterminée"


CONTRASTS = {                         # nom : (politique testée, référence)
    "dJ_R": ("causale", "r_fixe"),
    "dJ_M": ("causale", "m_fixe"),
    "dJ_tot": ("causale", "max(statique,directe)"),
    "clair_vs_R": ("clairvoyante", "r_fixe"),
}


def _campaign(grid, seeds, epochs, out_dir, base_cfg: N1Config | None = None) -> dict:
    """
    Boucle de campagne (V1.3.3). Appelée par run_grid() avec les paramètres
    préenregistrés ; les tests ne l'appellent qu'en miniature, sans publier de valeur.
    Écrit : runs.csv (une ligne par configuration × graine × politique),
    contrasts.csv, summary.json, SHA256SUMS.txt.
    """
    import csv
    import json
    os.makedirs(out_dir, exist_ok=True)
    base = base_cfg or N1Config()
    runs = []
    J = {}                                    # (k, politique) -> liste par graine
    for k, params in enumerate(grid):
        cfg = copy.deepcopy(base)
        for key, val in params.items():
            setattr(cfg, key, val)
        cfg.episode_epochs = epochs
        for s in seeds:
            res = simulate(cfg, s)
            for pol, m in res.items():
                runs.append(dict(config=k, w0_gs_m=cfg.w0_gs_m, c_r=cfg.c_r, seed=s, policy=pol, **m))
                J.setdefault((k, pol), []).append(m["J"])
    with open(os.path.join(out_dir, "runs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(runs[0].keys()))
        w.writeheader()
        w.writerows(runs)
    stats, sig, flags = {}, {}, {}
    for name, (pol, ref) in CONTRASTS.items():
        per = []
        for k in range(len(grid)):
            a = np.array(J[(k, pol)])
            b = (np.maximum(J[(k, "statique")], J[(k, "directe")]) if ref.startswith("max")
                 else np.array(J[(k, ref)]))
            per.append(contrast_stats(a, b))
        stats[name] = per
        sig[name] = significance_family(per)
        flags[name] = asymmetry_flags(per, sig[name])
    with open(os.path.join(out_dir, "contrasts.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["contrast", "config", "rel", "ci_low", "ci_high", "p_exact", "skew", "significant", "a_confirmer"])
        for name in CONTRASTS:
            for k, c in enumerate(stats[name]):
                w.writerow([name, k, c["rel"], c["ci"][0], c["ci"][1], c["p"], c["skew"],
                            bool(sig[name][k]), bool(flags[name][k])])
    issue = classify(sig["dJ_R"], sig["dJ_M"], sig["dJ_tot"], sig["clair_vs_R"])
    used = {"A": ["dJ_R", "dJ_tot"], "A partielle": ["dJ_R"], "A non exploitable": ["clair_vs_R"],
            "B": ["dJ_M"]}.get(issue, [])
    to_confirm = any(flags[n].any() for n in used)
    summary = dict(version="N1-01 V1.3.3", issue=issue, a_confirmer=to_confirm,
                   grid=[{k: float(v) for k, v in g.items()} for g in grid], seeds=list(seeds),
                   epochs=epochs, code_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                   v19_sha256=V19_SHA256, etsi_sha256=ETSI_PDF_SHA256,
                   significant={n: sig[n].tolist() for n in CONTRASTS},
                   a_confirmer_par_contraste={n: flags[n].tolist() for n in CONTRASTS})
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    with open(os.path.join(out_dir, "SHA256SUMS.txt"), "w") as f:
        for fn in ("runs.csv", "contrasts.csv", "summary.json"):
            f.write(hashlib.sha256(open(os.path.join(out_dir, fn), "rb").read()).hexdigest() + "  " + fn + "\n")
    return summary


def run_grid(authorization: str | None = None, out_dir: str = "n1_campagne_V1_3_3"):
    """Campagne unique préenregistrée : GRID (4) × SEEDS (20) × 1 000 époques."""
    if authorization != AUTH_TOKEN:
        raise PermissionError("Exécution N1-01 V1.3.3 non autorisée (§10).")
    if os.path.exists(os.path.join(out_dir, "summary.json")):
        raise RuntimeError("Campagne déjà exécutée dans ce répertoire : une seule campagne autorisée.")
    return _campaign(GRID, SEEDS, N1Config().episode_epochs, out_dir)
