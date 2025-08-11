# core.py
import numpy as np
import math
from scipy.linalg import eigvals
from scipy.integrate import solve_ivp
from .rates import k_marcus, k_eyring, k_tas_from_tau, kB_eV

class Reaction:
    def __init__(self, name, stoich, rtype, params, global_T):
        self.name = name
        self.stoich = stoich
        self.rtype = rtype.lower()
        self.params = params
        self.T = float(params.get('T', global_T))
        self.k_fwd = None
        self.k_rev = None
        self.reversible = bool(params.get('reversible', False))
        self._compute_rates()

    def _compute_rates(self):
        if self.rtype == 'marcus':
            V = float(self.params['V'])
            lam = float(self.params['lambda'])
            dG = float(self.params['dG'])
            kf = k_marcus(V, lam, dG, self.T)
            self.k_fwd = kf
            if self.reversible:
                K = math.exp(-dG / (kB_eV * self.T))
                self.k_rev = float(kf / K)

        elif self.rtype == 'eyring':
            dG_act = float(self.params['dG_act'])
            kf = k_eyring(dG_act, self.T)
            self.k_fwd = kf
            if self.reversible and 'dG' in self.params:
                dG = float(self.params['dG'])
                K = math.exp(-dG / (kB_eV * self.T))
                self.k_rev = float(kf / K)

        elif self.rtype == 'tas':
            tau = float(self.params['tau'])
            kf = k_tas_from_tau(tau)
            self.k_fwd = kf
            if self.reversible and 'dG' in self.params:
                dG = float(self.params['dG'])
                K = math.exp(-dG / (kB_eV * self.T))
                self.k_rev = float(kf / K)

        elif self.rtype == 'custom':
            if 'k_forward' in self.params:
                self.k_fwd = float(self.params['k_forward'])
            if 'k_reverse' in self.params:
                self.k_rev = float(self.params['k_reverse'])

    def rate(self, concentrations):
        rate_f = 1.0
        for sp, v in self.stoich.items():
            if v < 0:
                rate_f *= concentrations.get(sp, 0.0) ** abs(v)
        r_f = (self.k_fwd or 0.0) * rate_f

        r_r = 0.0
        if self.k_rev is not None:
            rate_r = 1.0
            for sp, v in self.stoich.items():
                if v > 0:
                    rate_r *= concentrations.get(sp, 0.0) ** v
            r_r = self.k_rev * rate_r
        return r_f, r_r

class MicrokineticModel:
    def __init__(self, species, reactions):
        self.species = species
        self.reactions = reactions
        self.idx = {s: i for i, s in enumerate(species)}
        self.S = self._build_stoich_matrix()

    def _build_stoich_matrix(self):
        S = np.zeros((len(self.species), len(self.reactions)))
        for j, r in enumerate(self.reactions):
            for sp, v in r.stoich.items():
                S[self.idx[sp], j] = v
        return S

    def rates_vector(self, y):
        conc = {s: float(y[self.idx[s]]) for s in self.species}
        rates = np.zeros(len(self.reactions))
        for j, r in enumerate(self.reactions):
            rf, rr = r.rate(conc)
            rates[j] = rf - rr
        return rates

    def rhs(self, t, y):
        return self.S.dot(self.rates_vector(y))

    def jacobian_approx(self, y0, eps=1e-8):
        n = len(y0)
        J = np.zeros((n, n))
        f0 = self.rhs(0.0, y0)
        for i in range(n):
            y_pert = y0.copy()
            dy = eps * max(1.0, abs(y0[i]))
            y_pert[i] += dy
            f1 = self.rhs(0.0, y_pert)
            J[:, i] = (f1 - f0) / dy
        return J

def detect_stiffness(model, y0):
    J = model.jacobian_approx(y0)
    eigs = eigvals(J)
    reals = np.abs(np.real(eigs))
    reals_nonzero = reals[reals > 1e-16]
    if reals_nonzero.size == 0:
        return 0.0, 0.0
    lam_max = reals_nonzero.max()
    lam_min = reals_nonzero.min()
    stiffness_ratio = float(lam_max / lam_min) if lam_min > 0 else float('inf')
    return float(lam_max), stiffness_ratio

def solve_model(model, y0, t_span, npoints):
    t_eval = np.linspace(t_span[0], t_span[1], npoints)
    lam_max, stiff_ratio = detect_stiffness(model, y0)
    method = 'BDF' if stiff_ratio > 1e3 else 'RK45'
    sol = solve_ivp(model.rhs, t_span, y0, method=method, t_eval=t_eval, atol=1e-12, rtol=1e-8)
    if not sol.success and method != 'BDF':
        sol = solve_ivp(model.rhs, t_span, y0, method='BDF', t_eval=t_eval)
    return sol

