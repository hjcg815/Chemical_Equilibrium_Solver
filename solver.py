# solver.py
import numpy as np
from thermo_data import thermo_data
from reactions import reactions  # your 20 reactions

R = 8.314  # J/mol·K
R_kJ = R / 1000  # kJ/mol·K

def calculate_reaction_thermo(reaction, T):
    """Calculate ΔH°, ΔS°, ΔG°, K for a reaction at temperature T"""
    ΔH = 0.0
    ΔS = 0.0

    for species, coeff in reaction['stoichiometry'].items():
        if species not in thermo_data:
            raise ValueError(f"{species} not in thermo_data.py")
        ΔH += coeff * thermo_data[species]['ΔHf']
        ΔS += coeff * thermo_data[species]['S'] / 1000  # J → kJ

    ΔG = ΔH - T * ΔS
    K = np.exp(-ΔG / (R_kJ * T))
    return ΔH, ΔS, ΔG, K

def solve_equilibrium(reaction, n0, T, P):
    """Solve for extent of reaction ξ and equilibrium mole fractions"""
    species_list = list(reaction['stoichiometry'].keys())
    nu = np.array([reaction['stoichiometry'][s] for s in species_list], dtype=float)
    n0_arr = np.array([n0.get(s, 0.0) for s in species_list], dtype=float)

    def residual(xi):
        n_eq = n0_arr + nu * xi
        if np.any(n_eq < 0):
            return 1e20
        Q = 1.0
        for i, ni in enumerate(n_eq):
            if nu[i] > 0:
                Q *= ni**nu[i] if ni > 0 else 1e-20
            elif nu[i] < 0:
                Q /= ni**(-nu[i]) if ni > 0 else 1e-20
        _, _, _, K = calculate_reaction_thermo(reaction, T)
        return K - Q

    xi_min = 0.0
    xi_max = 1e6
    tol = 1e-6
    max_iter = 1000
    xi = None
    for _ in range(max_iter):
        xi_mid = (xi_min + xi_max) / 2
        res = residual(xi_mid)
        if abs(res) < tol:
            xi = xi_mid
            break
        elif res > 0:
            xi_min = xi_mid
        else:
            xi_max = xi_mid
    if xi is None:
        xi = xi_mid

    n_eq = n0_arr + nu * xi
    n_total = np.sum(n_eq)
    y_eq = n_eq / n_total if n_total > 0 else np.zeros_like(n_eq)
    ΔH, ΔS, ΔG, K = calculate_reaction_thermo(reaction, T)

    results = {
        "ξ_eq": xi,
        "n_eq": dict(zip(species_list, n_eq)),
        "y_eq": dict(zip(species_list, y_eq)),
        "ΔH": ΔH,
        "ΔS": ΔS,
        "ΔG": ΔG,
        "K": K
    }
    return results
