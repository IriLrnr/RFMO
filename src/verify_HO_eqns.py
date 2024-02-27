# verify that all of my equations for the collective-harvesting optimum are giving the same result

import numpy as np
from scipy.optimize import root
from scipy.optimize import brentq
from my_functions import left_root_NE_fnc, dUidH_fixed_gamma_fnc

# ---
# parameters
# ---

# number of values to compare to verify the maxima
ngrid = 101

# use the default parameter values from the manuscript
pars = {
    "c": 0.05,  # cost of production per unit effort
    "wV": [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95],  # conservation weightings
}


# ---
# extract parameters
# ---

c = pars["c"]
wV = pars["wV"]
n = len(wV)

# ---
# find the historical allocations gamma
# ---

# find the historical harvesting effort
res = root(lambda hV: left_root_NE_fnc(pars, hV), [0.01] * n)
hsV = list(res.x)

# each player's ratio gamma_i^* = h_i / H is kept constant
Hs = sum(hsV)
gammaV = [hs / Hs for hs in hsV]

# ---
# method 1: find the root of d U_i / d H for all i
# ---

# this is the simplest method, which I'll use to verify the others below

HhatV = [None] * n
for i in range(n):
    HhatV[i] = brentq(lambda H: dUidH_fixed_gamma_fnc(pars, gammaV, i, H), 1e-6, Hs)

# ---
# verify that method 1 is indeed the maxima
# ---

print("")
print("== checking simplest method using derivative ==")
print("")

# list of U_i(\hat{H}) values
U_isV = [
    (1 - Hhat) ** wi * (gammai * Hhat * (1 - Hhat - c)) ** (1 - wi)
    for wi, gammai, Hhat in zip(wV, gammaV, HhatV)
]

for i in range(n):
    # vary H between its maximum and minimum
    HV = np.linspace(0, 1 - c, ngrid)

    # get the list of U_i(H) values
    U_iV = [(1 - H) ** wV[i] * (gammaV[i] * H * (1 - H - c)) ** (1 - wV[i]) for H in HV]

    # verify that U_i(\hat{H}) >= U_i(H) for all H
    U_is = U_isV[i]
    check_list = [U_is >= U_i for U_i in U_iV]
    print(f"U_{i}(Ĥ_{i}) >= U_{i}(H) for all H is {all(check_list)}")

# ---
# verify that Hhat above matches the left root of the quadratic equation
# ---

# find \hat{H} as the left root
HhatV_left_root = [
    1
    / 2
    * (
        c
        + 2 * wV[i]
        + np.sqrt(4 * c * wV[i] ** 2 + c**2 - 8 * c * wV[i] + 2 * c + 1)
        - 3
    )
    / (wV[i] - 2)
    for i in range(n)
]

print("")
print("== checking left-root method ==")
print("")

# print checking of matches
for i in range(n):
    diff = HhatV[i] - HhatV_left_root[i]
    within_tol = -1e-10 < diff < 1e-10
    print(f"Ĥ found for player {i} is within tolerance: {within_tol}")
