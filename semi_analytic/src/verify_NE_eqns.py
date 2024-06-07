# verify that all of my equations for the NE are giving the same result

import numpy as np
from scipy.optimize import root
from my_functions import dUdhiV_unconstrained_fnc, UV_fnc

# ---
# parameters
# ---

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
# method 1: find the root of d U_i / d h_i for all i
# ---

# this is the simplest method, which I'll use to verify the others below

res = root(lambda hV: dUdhiV_unconstrained_fnc(pars, hV), [0.01] * n)
hsV = list(res.x)


# ---
# verify that method 1 is indeed the maxima
# ---

print("")
print("== checking simplest method using derivative ==")
print("")

# list of U_i(h_i^*, H_{-i}^*) values
U_isV = UV_fnc(pars, hsV)

for i in range(n):
    # H_{-i}^* of other players
    H_i = sum(hs for j, hs in enumerate(hsV) if i != j)

    # vary h_i between its maxima and minima
    h_iV = np.linspace(0, 1 - c - H_i, 100)

    # get the list of U_i(h_i, H_{-i}^*) values
    U_iV = list()
    for h_i in h_iV:
        # update hV
        hV = [hs for hs in hsV]
        hV[i] = h_i

        # calculate new U_i
        U_i = UV_fnc(pars, hV)[i]
        U_iV.append(U_i)

    # verify that U_i(h_i^*, H_{-i}^*) >= U_i(h_i, H_{-i}^*) for all h_i
    U_is = U_isV[i]
    check_list = [U_is >= U_i for U_i in U_iV]
    print(f"U_{i}(h_{i}*) >= U_{i}(h_{i}) for all h_{i} is {all(check_list)}")


# ---
# method 2: find the root using my quadratic equation
# ---

print("")
print("")
print("== checking the quadratic gives the same answer ==")
print("")


# define the function for the quadratic
def Q_NE_fnc(pars, hV):
    # extract parameters
    c = pars["c"]
    wV = pars["wV"]
    n = len(wV)

    Q_NE = [None] * n
    for i in range(n):
        # parameters for player i
        w_i = wV[i]
        h_i = hV[i]
        H_i = sum(hV) - h_i  # H_{-i} of other players

        Q_NE[i] = (
            h_i**2 * (2 - w_i)
            + h_i * (H_i * (3 - 2 * w_i) + c + 2 * w_i - 3)
            + (w_i - 1) * (c - 1 - H_i * (H_i + c - 2))
        )

    return Q_NE


# find the root of the quadratic
res = root(lambda hV: Q_NE_fnc(pars, hV), [0.01] * n)
hsV_quad = list(res.x)

# verify that the root found using the quadratic matches the root found using the full equation
for i in range(n):
    diff = hsV[i] - hsV_quad[i]
    within_tol = -1e-10 < diff < 1e-10
    print(f"equilibrium found for {i} are within tolerance: {within_tol}")


# ---
# verify the left-hand root
# ---

left_root_NE = [None] * n
for i in range(n):
    # parameters for player i
    w_i = wV[i]
    h_i = hsV[i]
    H_i = sum(hsV) - h_i  # H_{-i} of other players

    # coefficients of the Q_NE equation
    A = 2 - w_i
    B = H_i * (3 - 2 * w_i) + c + 2 * w_i - 3
    C = (w_i - 1) * (c - 1 - H_i * (H_i + c - 2))

    # quadratic equation for left root
    left_root_NE[i] = (1 / (2 * A)) * (-B - np.sqrt(B**2 - 4 * A * C))

# ---
# method 4: failsafe method using left-root equation
# ---

print("")
print("")
print("== check left-root method ==")
print("")


# the most failsafe method uses the left root equation
def left_root_NE_fnc(pars, hV):
    # extract parameters
    c = pars["c"]
    wV = pars["wV"]
    n = len(wV)

    eq0 = [None] * n
    for i in range(n):
        # parameters for player i
        w_i = wV[i]
        h_i = hV[i]
        H_i = sum(hV) - h_i  # H_{-i} of other players

        # coefficients of the Q_NE equation
        A = 2 - w_i
        B = H_i * (3 - 2 * w_i) + c + 2 * w_i - 3
        C = (w_i - 1) * (c - 1 - H_i * (H_i + c - 2))

        # quadratic equation for left root - h_i
        eq0[i] = (1 / (2 * A)) * (-B - (B**2 - 4 * A * C) ** (1 / 2)) - h_i

    return eq0


# find the roots simultaneously
res = root(lambda hV: left_root_NE_fnc(pars, hV), [0.01] * n)
hsV_left_root = list(res.x)

# check within tolerance
for i in range(n):
    diff = hsV[i] - hsV_left_root[i]
    within_tol = -1e-10 < diff < 1e-10
    print(f"equilibrium found for {i} are within tolerance: {within_tol}")
