# explore how an increase or decrease in H might affect each player

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from my_functions import left_root_NE_fnc, UV_fnc, Hhat_fnc

# ---
# parameters
# ---

# number of grid points for plotting
ngrid = 50

# parameters for the scenario
pars = {
    "wV": [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95],  # conservation weightings
    "c": 0.05,  # cost of production per unit effort
}

# where the result directory is
dir_results = "../results/"

# ---
# extract parameters
# ---

c = pars["c"]
wV = pars["wV"]
n = len(wV)

# ---
# historically determined harvest allocations
# ---

# each player's ratio gamma_i = h_i^0 / H^0 is kept constant

# open access
res = root(lambda hV: left_root_NE_fnc(pars, hV), [0.01] * n)
hsV = list(res.x)
Hs = sum(hsV)
gamma_OAV = [hs / Hs for hs in hsV]

# exploratory
gamma_EXV = [1 / n for i in range(n)]

# ---
# players' collective harvesting optima
# ---

HhatV = [Hhat_fnc(pars, i) for i in range(n)]

# ---
# how each player's utility changes as H varies
# ---

# create a grid with more points near the ends
H_lo = 0
H_hi = 1 - c
dgrid = 0.1
HV = (
    list(np.linspace(H_lo, H_lo + dgrid, ngrid))
    + list(np.linspace(H_lo + dgrid, H_hi - dgrid, ngrid))
    + list(np.linspace(H_hi - dgrid, H_hi, ngrid))
)

# we'll do this for both scenarios, open access and explorator
gammaM = [gamma_OAV, gamma_EXV]
fnameV = ["H_deviation_open_access.pdf", "H_deviation_exploratory.pdf"]
H_histV = [Hs, n * 1e-3]
H_hist_labelV = [r"$H^*$", r"$H^{(0)}$"]
titleV = ["historically open access", "historically exploratory"]

for gammaV, fname, H_hist, H_hist_label, title in zip(
    gammaM, fnameV, H_histV, H_hist_labelV, titleV
):
    UM = [None] * ngrid * 3
    for idx, H_new in enumerate(HV):
        # the new hi's will be h_i = gamma_i * H_new
        hV = [gamma * H_new for gamma in gammaV]

        # vector of utilities (u_1, ..., u_n)
        UV = UV_fnc(pars, hV)

        # store
        UM[idx] = UV

    # plot

    plt.figure(figsize=(5, 3.8))
    plt.axvline(H_hist, ls="dashed", color="black", label=H_hist_label)

    UNM = np.array(UM)
    for foc in range(n):
        # how focal player's utility varies with H
        UV = UNM[:, foc]

        # plot how utility varies with collective-harvesting effort
        p = plt.plot(HV, UV, label=r"$i = " + str(foc + 1) + r"$")

        # plot a dot at focal player's optimum
        colour = p[0].get_color()
        Hhat = HhatV[foc]
        wi = wV[foc]
        hi = gammaV[foc] * Hhat
        U_Hhat = (1 - Hhat) ** wi * (hi * (1 - Hhat - c)) ** (1 - wi)
        plt.scatter([Hhat], [U_Hhat], color=colour)

    plt.title(title, fontsize="x-large")
    plt.legend(loc="best")
    plt.xlabel(r"collective-harvesting effort, $H$", fontsize="x-large")
    plt.ylabel(r"utility to focal player, $U_i(H)$", fontsize="x-large")
    plt.ylim((0, 0.8))  # so they can be compared easily
    plt.xlim((-0.01, 1.01))
    plt.tight_layout()
    plt.savefig(dir_results + fname)
    plt.close("all")
