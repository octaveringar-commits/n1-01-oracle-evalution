"""
fso_env.py — Environnement cross-layer Ground–LEO–GEO FSO (V19, étape 1, révision 2)
=====================================================================================

Décision à chaque époque Δt : a_t = (R_t, S_t, M_t)
    R_t : indice du relais LEO      (0..N-1)
    S_t : mode de relayage          (0 = AF, 1 = DF) — supprimé si s_fixed_df=True
    M_t : niveau MCS                (0..9  ↔ MCS 1..10)

Règles de provenance
--------------------
* Toutes les constantes sont déclarées dans Config, avec unité et source.
  « T4 » = Table 4 du manuscrit ; « H-xx » = hypothèse absente du manuscrit,
  soumise à validation de l'auteur AVANT tout entraînement (règle F).
* Aucun paramètre n'est ajusté pour retrouver des chiffres publiés antérieurement.
  Les valeurs par défaut des hypothèses H-xx sont celles déclarées en révision 1,
  ou des valeurs de fiche technique ; aucune n'a été choisie d'après le bilan.
* Dépendances : numpy, scipy. Gymnasium n'est requis que pour make_gym_env().

Révision 2 — corrections par rapport à la révision 1
----------------------------------------------------
C1  Détection cohérente (arbitrage D4) : SNR linéaire en puissance reçue,
    bruit de grenaille de l'oscillateur local + thermique [Ki16, Éq. 54 et 72].
    L'IM/DD reste disponible (detection="imdd") comme référence historique.
C2  Bande de bruit B déduite du débit symbole : B = R_s/2 (filtre adapté,
    bande 1/(2T) [Ki16]). R_s = 2 GBd reproduit la valeur B = 1 GHz déclarée
    en révision 1 (H-RS, en attente de validation).
C3  Relais AF : γ1 propre au relais AF (bruit ASE de l'EDFA, n_sp déduit du
    facteur de bruit) et loi de gain déclarée par κ (dérivation E5) :
        γ_AF = γ1,AF γ2 / (γ1,AF + γ2 + κ),  κ = 1 ⇔ Éq. (4) du manuscrit.
    Relais DF : γ1,DF = SNR du récepteur cohérent embarqué ; γ_DF = min(γ1,DF, γ2).
C4  Choix AF/DF : option s_fixed_df (S = DF imposé, action (R, M)) ou coût
    beta4_df ; les deux restent des décisions de l'auteur (D1, règle F).
C5  Éq. (14) appliquée strictement par défaut (reward_mode="indicator") ;
    la variante par fraction de sous-créneaux (H-BLK) est une option déclarée.
C6  Saut 1 (montant) : variance de scintillation onde sphérique émise au sol,
    intégrale exacte [z(L−z)/L]^{5/6} ; la formule de la rév. 1 (onde plane
    descendante) en est l'approximation (L−z)/L ≈ 1, qui surestime σ² d'≈ 1 %.

Références
----------
[AP05] L. C. Andrews, R. L. Phillips, Laser Beam Propagation through Random
       Media, 2nd ed., SPIE, 2005.
[FH07] A. A. Farid, S. Hranilovic, JLT 25(7):1702–1710, 2007.
[Ki16] K. Kikuchi, "Fundamentals of coherent optical fiber communications",
       JLT 34(1):157–179, 2016, DOI 10.1109/JLT.2015.2463719.
[Es10] R.-J. Essiambre et al., "Capacity limits of optical fiber networks",
       JLT 28(4):662–701, 2010, DOI 10.1109/JLT.2009.2039464.
[Ba12] K. Bayaki et al., IEEE TCOM 60(12):3797–3807, 2012 (loi de gain AF).

Limites connues (non corrigées, à déclarer dans le manuscrit)
-------------------------------------------------------------
L1  Paramètres α, β Gamma-Gamma : formules onde plane [AP05, chap. 9] appliquées
    à la variance montante ; approximation.
L2  Errance du faisceau montant incluse seulement via l'étalement long terme
    θ_eff (H-SPREAD), pas comme processus aléatoire.
L3  Normalisation γ de l'Éq. (9) : surrogate du manuscrit, non dérivé ici.
L4  γ1,AF dérivé pour un récepteur hétérodyne (E5, B = 2B_e) ; le cas
    homodyne + AF n'est pas dérivé [ND].
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfcinv, erfc

# ─────────────────────────────────────────────────────────────────────
# Constantes physiques (CODATA 2018)
# ─────────────────────────────────────────────────────────────────────
Q_E = 1.602176634e-19        # C
H_PLANCK = 6.62607015e-34    # J·s
C_LIGHT = 299_792_458.0      # m/s
R_EARTH_M = 6_371.0e3        # m (rayon moyen)
MU_EARTH = 3.986004418e14    # m^3/s^2

DETECTIONS = ("coherent_heterodyne", "coherent_homodyne", "imdd")


@dataclass
class Config:
    # ── Optique (Table 4) ──────────────────────────────────────────
    lambda_m: float = 1550e-9          # T4
    pt_gs_w: float = 1.0               # T4 : 30 dBm
    pt_leo_w: float = 1.0              # H-P1 : T4 ne distingue pas sol/LEO
    d_tx_m: float = 0.20               # T4 : ouverture d'émission (sol et LEO)
    d_rx_m: float = 0.30               # T4 : ouverture de réception (LEO et GEO)
    resp_a_per_w: float = 0.9          # T4 : responsivité (η_q = R·hν/q ≈ 0,72)
    n0_a2_per_hz: float = 1e-22        # T4 « 1e-22 W/Hz » ; H-N0 : DSP de bruit
                                       # thermique en courant (A²/Hz)
    eta_tx: float = 0.8                # H-OPT : rendement optique émission
    eta_rx: float = 0.8                # H-OPT : rendement optique réception

    # ── Détection (C1, C2) ─────────────────────────────────────────
    detection: str = "coherent_heterodyne"  # D4 ; voir DETECTIONS
    p_lo_w: float = 10e-3              # H-LO : oscillateur local (régime limité
                                       # par la grenaille de l'OL si 2qR·P_LO ≫ N0)
    symbol_rate_baud: float = 2.0e9    # H-RS : B = R_s/2 = 1 GHz (valeur de la rév. 1)
    include_signal_shot: bool = True   # grenaille du signal (négligeable en cohérent)

    # ── Relais AF (C3) ─────────────────────────────────────────────
    af_nf_db: float = 4.5              # H-NF : facteur de bruit EDFA spatial
                                       # (fiche Agiltron 4,5 dB ; autre source 6 dB)
    af_gain_db: float = 30.0           # H-G : gain EDFA (n'intervient que via (G−1)/G)
    af_kappa: float = 1.0              # H-KAPPA : 1 = Éq. (4) ; APC puissance totale
                                       # sur B_o ⇒ κ = B_o/B (E5) ; à décider
    af_noise_bw_factor: float = 2.0    # E5 : γ1,AF = P_r1 / (N_eq · 2B), hétérodyne

    # ── Géométrie ──────────────────────────────────────────────────
    h_leo_m: float = 550e3             # T4
    h_geo_m: float = 35_786e3          # T4
    gs_lat_deg: float = 12.13          # H-GS : N'Djamena (FTTJ)
    gs_alt_m: float = 295.0            # H-GS
    zenith_max_deg: float = 85.0       # T4 : plage 0–85°
    n_relays: int = 10                 # T4 : N ∈ [5, 20]
    dt_s: float = 0.1                  # T4 : époque de décision
    episode_steps: int = 1000          # 100 s par épisode

    # ── Pointage ───────────────────────────────────────────────────
    sigma_pe_range_rad: tuple = (2e-6, 8e-6)   # T4 : jitter sol→LEO, tiré par épisode
    sigma_pe_isl_rad: float = 1e-6             # H-PE2 : jitter LEO→GEO

    # ── Atmosphère ─────────────────────────────────────────────────
    cn2_ground: float = 1.7e-14        # HV5/7 : paramètre A (m^-2/3)
    hv_wind_rms_mps: float = 21.0      # HV5/7
    tau_zenith: float = 0.10           # H-ATM : épaisseur optique zénithale, ciel clair
    h_turb_top_m: float = 30e3         # borne d'intégration du profil HV

    # ── Statistique intra-époque ───────────────────────────────────
    n_subslots: int = 64               # H-BLK : sous-créneaux de cohérence par époque

    # ── ACM et récompense (Éq. 9 et 14) ────────────────────────────
    target_ber: float = 1e-9
    eta_m: float = 0.80                # T3 révisée
    rm: tuple = (0.50, 0.67, 0.80, 1.00, 1.20, 2.00, 3.00, 3.20, 4.00, 4.50)  # Table 2
    beta1: float = 10.0                # Éq. 14
    beta2: float = 0.1                 # Éq. 14 (MCS)
    beta3: float = 0.5                 # Éq. 14 (relais)
    reward_mode: str = "indicator"     # C5 : "indicator" = Éq. 14 ; "fraction" = H-BLK
    s_fixed_df: bool = False           # C4 : True ⇒ S = DF imposé (décision auteur)
    beta4_df: float = 0.0              # C4 : coût DF (latence/énergie) ; 0 = absent

    seed: int = 0

    @property
    def bw_hz(self) -> float:
        """Bande de bruit électrique B = R_s/2 (C2)."""
        return self.symbol_rate_baud / 2.0

    @property
    def photon_energy_j(self) -> float:
        return H_PLANCK * C_LIGHT / self.lambda_m


CONFIG = Config()


def validate_config(cfg: Config) -> None:
    if cfg.detection not in DETECTIONS:
        raise ValueError(f"detection doit appartenir à {DETECTIONS}")
    if cfg.reward_mode not in ("indicator", "fraction"):
        raise ValueError("reward_mode doit valoir 'indicator' ou 'fraction'")
    if cfg.af_kappa < 1.0:
        raise ValueError("κ = B_o/B ≥ 1 (E5)")
    if 10 ** (cfg.af_nf_db / 10) < 1.0:
        raise ValueError("facteur de bruit < 0 dB impossible")


# ─────────────────────────────────────────────────────────────────────
# ACM : Éq. 9 et Table 3 révisée
# ─────────────────────────────────────────────────────────────────────
def gamma_th_linear(cfg: Config) -> np.ndarray:
    """Seuils γ_th,m tels que BER_m(γ_th,m) = target_ber exactement (Éq. 9)."""
    rm = np.asarray(cfg.rm)
    x = erfcinv(2.0 * cfg.target_ber)
    return x ** 2 * (2.0 ** rm - 1.0) / cfg.eta_m


def ber_eq9(gamma_lin, rm, eta_m: float):
    """Éq. 9 : BER_m(γ) = ½ erfc( sqrt(γ η_m / (2^R_m − 1)) )."""
    arg = np.sqrt(np.maximum(gamma_lin, 0.0) * eta_m / (2.0 ** np.asarray(rm) - 1.0))
    return 0.5 * erfc(arg)


# ─────────────────────────────────────────────────────────────────────
# Turbulence : HV5/7, Rytov, Fried, Gamma-Gamma
# ─────────────────────────────────────────────────────────────────────
def cn2_hv(h_m: float, cfg: Config) -> float:
    """Profil Hufnagel–Valley généralisé [AP05, chap. 12], h au-dessus du sol (m)."""
    v = cfg.hv_wind_rms_mps
    return (0.00594 * (v / 27.0) ** 2 * (1e-5 * h_m) ** 10 * np.exp(-h_m / 1000.0)
            + 2.7e-16 * np.exp(-h_m / 1500.0)
            + cfg.cn2_ground * np.exp(-h_m / 100.0))


def slant_range_leo(zenith_rad: float, cfg: Config) -> float:
    """Distance oblique station → LEO (Terre sphérique) pour l'angle zénithal θ."""
    r0 = R_EARTH_M + cfg.gs_alt_m
    rs = R_EARTH_M + cfg.h_leo_m
    c = np.cos(zenith_rad)
    return float(np.sqrt((r0 * c) ** 2 + rs ** 2 - r0 ** 2) - r0 * c)


def rytov_variance(zenith_rad: float, cfg: Config) -> float:
    """
    Variance de scintillation, saut 1 montant (sol → LEO), onde sphérique émise
    au sol [AP05, chap. 5 et 12] (C6) :
        σ² = 2,25 k^{7/6} ∫_0^L C_n²(z) [z (L − z)/L]^{5/6} dz,   z = h·sec θ (H-FLAT)
    L = distance oblique sur Terre sphérique. Forme exacte utilisée par
    l'environnement ; l'approximation (L − z)/L ≈ 1 (formule onde plane
    descendante) surestime σ² d'environ 1 % : voir rytov_variance_flat_approx().
    """
    k = 2.0 * np.pi / cfg.lambda_m
    sec = 1.0 / np.cos(zenith_rad)
    L = slant_range_leo(zenith_rad, cfg)

    def f(h):
        z = h * sec
        return cn2_hv(h, cfg) * (z * (L - z) / L) ** (5.0 / 6.0) * sec

    integral, _ = quad(f, 0.0, cfg.h_turb_top_m, limit=200, points=[100, 1000, 10000])
    return 2.25 * k ** (7.0 / 6.0) * integral


def rytov_variance_flat_approx(zenith_rad: float, cfg: Config) -> float:
    """Approximation (L − z)/L ≈ 1 (formule de la rév. 1), conservée pour contrôle :
        σ² ≈ 2,25 k^{7/6} sec^{11/6}(θ) ∫ C_n²(h) h^{5/6} dh."""
    k = 2.0 * np.pi / cfg.lambda_m
    integral, _ = quad(lambda h: cn2_hv(h, cfg) * h ** (5.0 / 6.0),
                       0.0, cfg.h_turb_top_m, limit=200, points=[100, 1000, 10000])
    return 2.25 * k ** (7.0 / 6.0) * (1.0 / np.cos(zenith_rad)) ** (11.0 / 6.0) * integral


def fried_parameter(zenith_rad: float, cfg: Config) -> float:
    """
    r0 = [0,423 k² sec θ ∫ C_n²(h) w(h) dh]^{-3/5} [AP05, chap. 12].
    Montant, onde sphérique : w(h) = (1 − h/H_sat)^{5/3} ≈ 1 pour h ≤ 30 km ≪ 550 km.
    """
    k = 2.0 * np.pi / cfg.lambda_m
    integral, _ = quad(lambda h: cn2_hv(h, cfg) * (1.0 - h / cfg.h_leo_m) ** (5.0 / 3.0),
                       0.0, cfg.h_turb_top_m, limit=200, points=[100, 1000, 10000])
    return (0.423 * k ** 2 * integral / np.cos(zenith_rad)) ** (-3.0 / 5.0)


def gamma_gamma_params(sigma_r2: float) -> tuple[float, float]:
    """
    α, β Gamma-Gamma, onde plane, échelle interne nulle [AP05, chap. 9] :
      α = [exp(0,49 σ_R² / (1 + 1,11 σ_R^{12/5})^{7/6}) − 1]^{-1}
      β = [exp(0,51 σ_R² / (1 + 0,69 σ_R^{12/5})^{5/6}) − 1]^{-1}
    σ_R^{12/5} = (σ_R²)^{6/5}. Écrire sigma_r2**(12/5) serait faux.
    """
    s = max(sigma_r2, 1e-9)
    s125 = s ** (6.0 / 5.0)
    alpha = 1.0 / (np.exp(0.49 * s / (1.0 + 1.11 * s125) ** (7.0 / 6.0)) - 1.0)
    beta = 1.0 / (np.exp(0.51 * s / (1.0 + 0.69 * s125) ** (5.0 / 6.0)) - 1.0)
    return alpha, beta


def sample_gamma_gamma(alpha: float, beta: float, size, rng: np.random.Generator):
    """I_a = X·Y, X ~ Γ(α, 1/α), Y ~ Γ(β, 1/β) ⇒ E[I_a] = 1."""
    return rng.gamma(alpha, 1.0 / alpha, size) * rng.gamma(beta, 1.0 / beta, size)


# ─────────────────────────────────────────────────────────────────────
# Pertes de pointage [FH07]
# ─────────────────────────────────────────────────────────────────────
def pointing_params(w_m, a_m: float):
    """A0 = erf(v)², w_eq² = w² √π erf(v) / (2 v e^{−v²}), v = √π a / (√2 w)."""
    v = np.sqrt(np.pi) * a_m / (np.sqrt(2.0) * w_m)
    a0 = erf(v) ** 2
    weq2 = w_m ** 2 * np.sqrt(np.pi) * erf(v) / (2.0 * v * np.exp(-v ** 2))
    return a0, weq2


# ─────────────────────────────────────────────────────────────────────
# Géométrie : repère terrestre fixe (rotation de la Terre négligée sur 100 s, H-ROT)
# ─────────────────────────────────────────────────────────────────────
def gs_position(cfg: Config) -> np.ndarray:
    lat = np.radians(cfg.gs_lat_deg)
    r = R_EARTH_M + cfg.gs_alt_m
    return r * np.array([np.cos(lat), 0.0, np.sin(lat)])


def geo_position(cfg: Config) -> np.ndarray:
    """GEO à l'équateur, même longitude que la station (H-GEO)."""
    return (R_EARTH_M + cfg.h_geo_m) * np.array([1.0, 0.0, 0.0])


def zenith_angle(gs: np.ndarray, sat: np.ndarray) -> np.ndarray:
    up = gs / np.linalg.norm(gs)
    vec = sat - gs
    cosz = (vec @ up) / np.linalg.norm(vec, axis=-1)
    return np.arccos(np.clip(cosz, -1.0, 1.0))


def line_of_sight_clear(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    """Vrai si le segment p1→p2 n'intercepte pas la Terre (rayon R_EARTH)."""
    d = p2 - p1
    t = np.clip(-(p1 * d).sum(-1) / (d * d).sum(-1), 0.0, 1.0)
    closest = p1 + t[..., None] * d
    return np.linalg.norm(closest, axis=-1) > R_EARTH_M


def leo_at_zenith(zdeg: float, cfg: Config) -> np.ndarray:
    """Point LEO à l'angle zénithal z, dans le plan méridien, côté équateur."""
    gs = gs_position(cfg)
    up = gs / np.linalg.norm(gs)
    south = np.cross(np.array([0.0, 1.0, 0.0]), up)
    south /= np.linalg.norm(south)
    z = np.radians(zdeg)
    dvec = np.cos(z) * up - np.sin(z) * south
    r_leo = R_EARTH_M + cfg.h_leo_m
    b = 2 * gs @ dvec
    c = gs @ gs - r_leo ** 2
    return gs + (-b + np.sqrt(b * b - 4 * c)) / 2 * dvec


class Constellation:
    """
    N relais LEO sur orbites circulaires (h_leo), grands cercles d'orientation
    aléatoire (H-ORB). Position : p(t) = r [cos(ωt) u0 + sin(ωt) v0].
    Chaque relais est tiré visible à t=0 (zénith < 80°).
    """

    def __init__(self, cfg: Config, rng: np.random.Generator):
        self.cfg = cfg
        r = R_EARTH_M + cfg.h_leo_m
        self.r = r
        self.omega = np.sqrt(MU_EARTH / r ** 3)
        gs = gs_position(cfg)
        u0s, v0s = [], []
        while len(u0s) < cfg.n_relays:
            u = rng.normal(size=3)
            u /= np.linalg.norm(u)
            if np.degrees(zenith_angle(gs, r * u)) < 80.0:
                w = rng.normal(size=3)
                w -= (w @ u) * u
                w /= np.linalg.norm(w)
                u0s.append(u)
                v0s.append(w)
        self.u0 = np.array(u0s)
        self.v0 = np.array(v0s)

    def positions(self, t_s: float) -> np.ndarray:
        c, s = np.cos(self.omega * t_s), np.sin(self.omega * t_s)
        return self.r * (c * self.u0 + s * self.v0)


# ─────────────────────────────────────────────────────────────────────
# Bilan de liaison
# ─────────────────────────────────────────────────────────────────────
def slant_range(gs: np.ndarray, sats: np.ndarray) -> np.ndarray:
    return np.linalg.norm(sats - gs, axis=-1)


def receiver_snr(pr_w, cfg: Config):
    """
    SNR électrique du récepteur pour une puissance optique reçue P_r (C1).

    Cohérent, OL de puissance P_LO, bande de bruit B = R_s/2 [Ki16] :
      homodyne   : γ = 4 R² P_LO P_r / σ²
      hétérodyne : γ = 2 R² P_LO P_r / σ²
      σ² = 2 q R (P_LO + P_r) B + N0 B
    Limite grenaille de l'OL : γ_hom = 2 R P_r/(q B), γ_het = R P_r/(q B),
    soit 4 n_s et 2 n_s photons détectés par symbole [Ki16, Éq. 54 et 72].

    IM/DD (référence historique, non retenue — D4) :
      γ = (R P_r)² / (N0 B + 2 q R P_r B).
    """
    pr = np.asarray(pr_w, dtype=float)
    R, B, N0 = cfg.resp_a_per_w, cfg.bw_hz, cfg.n0_a2_per_hz
    if cfg.detection == "imdd":
        var = N0 * B + (2.0 * Q_E * R * pr * B if cfg.include_signal_shot else 0.0)
        return (R * pr) ** 2 / var
    p_shot = cfg.p_lo_w + (pr if cfg.include_signal_shot else 0.0)
    var = 2.0 * Q_E * R * p_shot * B + N0 * B
    mix = 4.0 if cfg.detection == "coherent_homodyne" else 2.0
    return mix * R ** 2 * cfg.p_lo_w * pr / var


def af_nsp(cfg: Config) -> float:
    """n_sp tel que NF = 2 n_sp (G−1)/G + 1/G [Es10, §IV ; Ki16, Éq. 85–86]."""
    nf = 10 ** (cfg.af_nf_db / 10)
    g = 10 ** (cfg.af_gain_db / 10)
    return (nf - 1.0 / g) * g / (2.0 * (g - 1.0))


def af_hop1_snr(pr_w, cfg: Config):
    """
    SNR du premier saut pour le relais AF tout-optique (C3, dérivation E5) :
      γ1,AF = P_r1 / (N_eq · 2B),  N_eq = n_sp hν (G−1)/G.
    En IM/DD (référence historique), γ1,AF = γ du récepteur (Éq. 4 d'origine).
    """
    pr = np.asarray(pr_w, dtype=float)
    if cfg.detection == "imdd":
        return receiver_snr(pr, cfg)
    g = 10 ** (cfg.af_gain_db / 10)
    n_eq = af_nsp(cfg) * cfg.photon_energy_j * (g - 1.0) / g
    return pr / (n_eq * cfg.af_noise_bw_factor * cfg.bw_hz)


def theta_diff(cfg: Config) -> float:
    """Divergence (demi-angle 1/e²) d'un faisceau gaussien w0 = D_tx/2."""
    return cfg.lambda_m / (np.pi * cfg.d_tx_m / 2.0)


def hop1_mean_power(zenith_rad, dist_m, cfg: Config, r0_m):
    """
    Puissance reçue moyenne (sans fading, sans erreur de pointage) sol → LEO.
    Divergence montante effective (étalement long terme, H-SPREAD) :
        θ_eff = sqrt(θ_diff² + (λ / r0)²)
    """
    th = np.sqrt(theta_diff(cfg) ** 2 + (cfg.lambda_m / r0_m) ** 2)
    w = th * dist_m
    a0, weq2 = pointing_params(w, cfg.d_rx_m / 2.0)
    latm = np.exp(-cfg.tau_zenith / np.cos(zenith_rad))
    pr = cfg.pt_gs_w * cfg.eta_tx * cfg.eta_rx * latm * a0
    return pr, weq2, latm, a0


def hop2_mean_power(dist_m, cfg: Config):
    """LEO → GEO, vide : divergence de diffraction seule."""
    w = theta_diff(cfg) * dist_m
    a0, weq2 = pointing_params(w, cfg.d_rx_m / 2.0)
    pr = cfg.pt_leo_w * cfg.eta_tx * cfg.eta_rx * a0
    return pr, weq2, a0


def combine_af_df(g1_af, g1_df, g2, kappa: float = 1.0):
    """
    AF (C3, E5) : γ1,AF γ2 / (γ1,AF + γ2 + κ) ; κ = 1 et γ1,AF = γ1,DF ⇔ Éq. 4.
    DF (Éq. 5)  : min(γ1,DF, γ2).
    """
    af = g1_af * g2 / (g1_af + g2 + kappa)
    df = np.minimum(g1_df, g2)
    return af, df


# ─────────────────────────────────────────────────────────────────────
# Environnement
# ─────────────────────────────────────────────────────────────────────
class CrossLayerFSOEnv:
    """
    API de type Gymnasium (reset/step) sans dépendance.

    Observation (H-OBS), par relais i, mesures de l'époque t−1 (CSI retardée, H-CSI) :
        [γ̄_1,DF,i (dB)/40, γ̄_2,i (dB)/40, θ_i/90°, visible_i]
      + global : [log10(C_n²_sol)+14, R_{t−1}/N, S_{t−1}, M_{t−1}/9]
    γ̄ = moyenne en dB sur les sous-créneaux de l'époque précédente.
    (γ1,AF est une fonction déterministe de γ1,DF : P_r1 commun ; non répété.)

    Récompense :
      reward_mode="indicator" (Éq. 14 stricte) :
        v = 1[ BER_moy > 1e-9 ],  BER_moy = moyenne de BER_m(γ_eff) sur l'époque
        r = R_m (1 − v) − β1 v − β2 1[M≠M⁻] − β3 1[R≠R⁻] − β4 1[S=DF]
      reward_mode="fraction" (H-BLK) : v = fraction de sous-créneaux en violation.
    Relais non visible ⇒ v = 1, débit nul.
    """

    def __init__(self, cfg: Config | None = None):
        self.cfg = copy.deepcopy(cfg or CONFIG)
        validate_config(self.cfg)
        self.rm = np.asarray(self.cfg.rm)
        self.gth = gamma_th_linear(self.cfg)
        self.n_obs = 4 * self.cfg.n_relays + 4
        n_modes = 1 if self.cfg.s_fixed_df else 2
        self.action_nvec = (self.cfg.n_relays, n_modes, len(self.rm))
        self.rng = np.random.default_rng(self.cfg.seed)
        self._turb_cache: dict = {}

    def _turb(self, zenith_rad: float):
        key = (round(np.degrees(zenith_rad) * 2) / 2, self.cfg.cn2_ground)
        if key not in self._turb_cache:
            z = np.radians(key[0])
            s2 = rytov_variance(z, self.cfg)
            r0 = fried_parameter(z, self.cfg)
            self._turb_cache[key] = (s2, *gamma_gamma_params(s2), r0)
        return self._turb_cache[key]

    def _draw_epoch(self, t_s: float):
        """Retourne zen, vis, P_r1 (N,K), P_r2 (N,K) en W."""
        cfg, rng, K = self.cfg, self.rng, self.cfg.n_subslots
        sats = self.const.positions(t_s)
        gs, geo = self.gs, self.geo
        zen = zenith_angle(gs, sats)
        vis = (np.degrees(zen) < cfg.zenith_max_deg) & line_of_sight_clear(gs, sats) \
            & line_of_sight_clear(sats, geo)
        N = cfg.n_relays
        p1 = np.zeros((N, K))
        p2 = np.zeros((N, K))
        d1 = slant_range(gs, sats)
        d2 = np.linalg.norm(geo - sats, axis=-1)
        for i in range(N):
            if not vis[i]:
                continue
            _, a, b, r0 = self._turb(zen[i])
            pr1, weq2, _, _ = hop1_mean_power(zen[i], d1[i], cfg, r0)
            ia = sample_gamma_gamma(a, b, K, rng)
            rad = rng.rayleigh(self.sigma_pe * d1[i], K)
            p1[i] = pr1 * ia * np.exp(-2.0 * rad ** 2 / weq2)
            pr2, weq2b, _ = hop2_mean_power(d2[i], cfg)
            rad2 = rng.rayleigh(cfg.sigma_pe_isl_rad * d2[i], K)
            p2[i] = pr2 * np.exp(-2.0 * rad2 ** 2 / weq2b)
        return zen, vis, p1, p2

    def snrs(self, p1, p2):
        """(γ1,AF, γ1,DF, γ2) à partir des puissances reçues."""
        return af_hop1_snr(p1, self.cfg), receiver_snr(p1, self.cfg), receiver_snr(p2, self.cfg)

    @staticmethod
    def _db_mean(g):
        return 10.0 * np.log10(np.maximum(g, 1e-6)).mean(axis=-1)

    def _obs(self):
        cfg = self.cfg
        _, g1d, g2 = self.snrs(self._p1, self._p2)
        per = np.stack([
            np.clip(self._db_mean(g1d) / 40.0, -2, 2) * self._vis,
            np.clip(self._db_mean(g2) / 40.0, -2, 2) * self._vis,
            np.degrees(self._zen) / 90.0,
            self._vis.astype(float),
        ], axis=1).ravel()
        glob = np.array([np.log10(cfg.cn2_ground) + 14.0,
                         self._prev[0] / cfg.n_relays, self._prev[1], self._prev[2] / 9.0])
        return np.concatenate([per, glob]).astype(np.float32)

    def reset(self, seed: int | None = None):
        if seed is not None:
            self.rng = np.random.default_rng(seed)
        cfg = self.cfg
        self.gs, self.geo = gs_position(cfg), geo_position(cfg)
        self.const = Constellation(cfg, self.rng)
        self.sigma_pe = self.rng.uniform(*cfg.sigma_pe_range_rad)
        self.t = 0
        self._zen, self._vis, self._p1, self._p2 = self._draw_epoch(0.0)
        _, g1d, _ = self.snrs(self._p1, self._p2)
        self._prev = (int(np.argmax(self._db_mean(g1d) * self._vis - 1e9 * (~self._vis))), 1, 0)
        return self._obs(), {}

    def _mode(self, mode: int) -> int:
        return 1 if self.cfg.s_fixed_df else int(mode)

    def evaluate(self, relay: int, mode: int, mcs: int, zen, vis, p1, p2):
        """Débit, violation v et γ_eff pour une action sur une époque tirée."""
        cfg = self.cfg
        if not vis[relay]:
            return 0.0, 1.0, np.zeros(cfg.n_subslots)
        g1a, g1d, g2 = self.snrs(p1[relay], p2[relay])
        af, df = combine_af_df(g1a, g1d, g2, cfg.af_kappa)
        geff = df if self._mode(mode) == 1 else af
        ber = ber_eq9(geff, self.rm[mcs], cfg.eta_m)
        if cfg.reward_mode == "indicator":
            v = float(ber.mean() > cfg.target_ber)
        else:
            v = float((ber > cfg.target_ber).mean())
        return float(self.rm[mcs] * (1.0 - v)), v, geff

    def step(self, action):
        cfg = self.cfg
        relay, mode, mcs = (int(x) for x in action)
        mode = self._mode(mode)
        self.t += 1
        zen, vis, p1, p2 = self._draw_epoch(self.t * cfg.dt_s)
        thr, v, geff = self.evaluate(relay, mode, mcs, zen, vis, p1, p2)
        sw_m = float(mcs != self._prev[2])
        sw_r = float(relay != self._prev[0])
        reward = (thr - cfg.beta1 * v - cfg.beta2 * sw_m - cfg.beta3 * sw_r
                  - cfg.beta4_df * float(mode == 1))
        self._zen, self._vis, self._p1, self._p2 = zen, vis, p1, p2
        self._prev = (relay, mode, mcs)
        info = {"se": thr, "violation": v, "switch_mcs": sw_m, "switch_relay": sw_r,
                "mode": mode, "zenith_deg": float(np.degrees(zen[relay])),
                "visible": bool(vis[relay]),
                "geff_db_mean": float(self._db_mean(geff[None])[0]) if vis[relay] else -60.0}
        truncated = self.t >= cfg.episode_steps
        return self._obs(), float(reward), False, truncated, info


def make_gym_env(cfg: Config | None = None):
    """Enveloppe Gymnasium (MultiDiscrete) pour l'entraînement PPO."""
    import gymnasium as gym
    from gymnasium import spaces

    class _GymEnv(gym.Env):
        metadata = {"render_modes": []}

        def __init__(self):
            super().__init__()
            self.core = CrossLayerFSOEnv(cfg)
            self.action_space = spaces.MultiDiscrete(list(self.core.action_nvec))
            self.observation_space = spaces.Box(-np.inf, np.inf, (self.core.n_obs,),
                                                np.float32)

        def reset(self, seed=None, options=None):
            super().reset(seed=seed)
            return self.core.reset(seed=seed)

        def step(self, action):
            return self.core.step(action)

    return _GymEnv()
