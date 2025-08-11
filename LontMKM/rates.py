# rates.py
import math

# Constants in eV units
kB_eV = 8.617333262145e-5
h_eVs = 4.135667696e-15
hbar_eVs = h_eVs / (2 * math.pi)

def k_marcus(V_eV, lam_eV, dG_eV, T):
    """Marcus theory rate constant."""
    prefactor = 2 * math.pi / hbar_eVs
    denom = math.sqrt(4 * math.pi * lam_eV * kB_eV * T)
    exponent = -((dG_eV + lam_eV) ** 2) / (4 * lam_eV * kB_eV * T)
    return float(prefactor * (V_eV ** 2) / denom * math.exp(exponent))

def k_eyring(dG_act_eV, T):
    """Eyring equation rate constant."""
    pre = (kB_eV * T) / h_eVs
    return float(pre * math.exp(-dG_act_eV / (kB_eV * T)))

def k_tas_from_tau(tau_s):
    """Rate constant from transient absorption spectroscopy (1/tau)."""
    if tau_s <= 0:
        raise ValueError("tau must be positive")
    return 1.0 / tau_s

