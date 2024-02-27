# solves and plots scenarios: historical x voting system

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from my_functions import left_root_NE_fnc, Hhat_fnc, UV_fnc

# ---
# parameters
# ---

# where the result directory is
dir_results = "../results/"

# parameter values for all scenarios
pars = {
    "c": 0.05,  # cost of production per unit effort
    "wV": [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95],  # conservation weightings
}

# historical scenarios
histV = ["open access", "exploratory"]

# voting scenarios, also orders voting systems for plotting
voting_systemV = ["consensus", "supermajority", "majority"]

# the proportion required for the vote to pass
vD = {"consensus": 1, "supermajority": 0.75, "majority": 0.5}

# colours for plotting
scenario2colour = {
    "exploratory": "#440154",
    "open access": "#440154",
    "consensus": "#3b528b",
    "supermajority": "#21908c",
    "majority": "#5dc963",
}

# ---
# extract parameters
# ---

wV = pars["wV"]
c = pars["c"]
n = len(wV)


# ---
# find the historical harvesting efforts and allocations
# ---

# open access historical harvesting effort is the Nash equilibrium
res = root(lambda hV: left_root_NE_fnc(pars, hV), [0.01] * n)
h_hist_OAV = list(res.x)
H_hist_OA = sum(h_hist_OAV)
gamma_OAV = [h / H_hist_OA for h in h_hist_OAV]
U_OAV = UV_fnc(pars, h_hist_OAV)
N_hist_OA = 1 - H_hist_OA

# exploratory historical harvesting effort is 1e-3 for every player
h_hist_EXV = [1e-3 for i in range(n)]

H_hist_EX = sum(h_hist_EXV)
gamma_EXV = [1 / n for i in range(n)]
U_EXV = UV_fnc(pars, h_hist_EXV)
N_hist_EX = 1 - H_hist_EX

# put into dictionaries
histD = {
    "open access": {
        "hV": h_hist_OAV,
        "H": H_hist_OA,
        "N": N_hist_OA,
        "gammaV": gamma_OAV,
        "UV": U_OAV,
    },
    "exploratory": {
        "hV": h_hist_EXV,
        "H": H_hist_EX,
        "N": N_hist_EX,
        "gammaV": gamma_EXV,
        "UV": U_EXV,
    },
}


# ---
# find the voting-dynamical endpoint
# ---

# the functions that define the two potential critical players
Klo_fnc = lambda v, n: int(np.ceil(v * n))  # noqa: E731
Khi_fnc = lambda v, n: n + 1 - int(np.ceil(n * v))  # noqa: E731

# take advantage of the fact that they are already ordered in order of increasing environmental concern
KloD = {vot_sys: Klo_fnc(vD[vot_sys], n) for vot_sys in voting_systemV}
KhiD = {vot_sys: Khi_fnc(vD[vot_sys], n) for vot_sys in voting_systemV}

# another version as python indices count from 0
pKloD = {vot_sys: Klo - 1 for vot_sys, Klo in KloD.items()}
pKhiD = {vot_sys: Khi - 1 for vot_sys, Khi in KhiD.items()}

# the Nash equilibrium is higher than any optimum, so for the open-access scenario, we use Khi
# the exporatory scenario starts at arbitrarily small harvesting, so we use Klo
pKD = {"open access": pKhiD, "exploratory": pKloD}

# find each player's collective-harvesting optimum
HhatV = [Hhat_fnc(pars, i) for i in range(n)]

# find the critical player's collective-harvesting optimum for each scenario
# this is the voting-dynamical endpoint
HvotD = {
    hist: {vot_sys: HhatV[pKD[hist][vot_sys]] for vot_sys in voting_systemV}
    for hist in histV
}

# ---
# find the characteristics of the voting-dynamical endpoints
# ---

# harvesting efforts
hvotD = {
    hist: {
        vot_sys: [gamma * HvotD[hist][vot_sys] for gamma in histD[hist]["gammaV"]]
        for vot_sys in voting_systemV
    }
    for hist in histV
}

# utilities def UV_fnc(pars, hV):
UvotD = {
    hist: {vot_sys: UV_fnc(pars, hvotD[hist][vot_sys]) for vot_sys in voting_systemV}
    for hist in histV
}

# abundances
NvotD = {
    hist: {vot_sys: 1 - HvotD[hist][vot_sys] for vot_sys in voting_systemV}
    for hist in histV
}

# ---
# plot each of the short plots
# ---

for hist in histV:
    # total harvest effort
    # ---

    ylabel = r"total harvesting effort $H$"

    HV = [histD[hist]["H"]] + [HvotD[hist][vot_sys] for vot_sys in voting_systemV]
    scenarioV = [hist] + voting_systemV
    colourV = [scenario2colour[scen] for scen in scenarioV]

    plt.figure(figsize=(5, 3.8))
    ax = plt.subplot(111)
    xticks = list(range(4))
    ax.bar(xticks, HV, color=colourV, align="center")
    ax.set_xticks(xticks)
    ax.set_xticklabels(scenarioV)
    ax.set_ylabel(ylabel, fontsize="x-large")
    ax.set_ylim((0, 1))
    plt.tight_layout()
    fname = dir_results + "scenario_" + hist.replace(" ", "_") + "_H.pdf"
    plt.savefig(fname)
    plt.close("all")

    # abundance
    # ---

    ylabel = r"abundance $n^*$"

    NV = [histD[hist]["N"]] + [NvotD[hist][vot_sys] for vot_sys in voting_systemV]
    scenarioV = [hist] + voting_systemV
    colourV = [scenario2colour[scen] for scen in scenarioV]

    plt.figure(figsize=(5, 3.8))
    ax = plt.subplot(111)
    xticks = list(range(4))
    ax.bar(xticks, NV, color=colourV, align="center")
    ax.set_xticks(xticks)
    ax.set_xticklabels(scenarioV)
    ax.set_ylabel(ylabel, fontsize="x-large")
    ax.set_ylim((0, 1))
    plt.tight_layout()
    fname = dir_results + "scenario_" + hist.replace(" ", "_") + "_N.pdf"
    plt.savefig(fname)
    plt.close("all")

# ---
# plot each of the long plots
# ---

# width and positions of bars
lv = len(voting_systemV) + 1
xdel = 1 / (lv + 1)
width = 1 / (lv + 2)
x = np.array(range(n))

for hist in histV:
    # harvesting efforts
    # ---

    plt.figure(figsize=(7, 3.8))
    ax = plt.subplot(111)
    ax.set_yticks([0, 0.1, 0.2])  # NOTE hardcoded for ease of comparison
    ax.set_ylim((0, 0.2))  # NOTE hardcoded for ease of comparison

    # store each bar chart element to add legends later
    rects = list()

    # plot historical scenario
    hV = histD[hist]["hV"]
    colour = scenario2colour[hist]
    rect = ax.bar(x, hV, width=width, color=colour, align="center")
    rects.append(rect[0])

    # plot each voting-dynamical endpoint
    for idx, vot_sys in enumerate(voting_systemV):
        hV = hvotD[hist][vot_sys]
        colour = scenario2colour[vot_sys]
        rect = ax.bar(
            x + (idx + 1) * xdel, hV, width=width, color=colour, align="center"
        )
        rects.append(rect[0])

    # deocrate the plot and save
    ylabel = r"harvesting efforts $h_i$"
    xticks_w = x + xdel * (lv - 1) / 2
    xticklabels_w = [
        r"$w_" + str(foc + 1) + " = " + str(wV[foc]) + "$" for foc in range(n)
    ]
    ax.set_ylabel(ylabel, fontsize="x-large")
    ax.set_xticks(xticks_w)
    ax.set_xticklabels(xticklabels_w)
    ax.legend(rects, [hist] + voting_systemV)
    plt.tight_layout()
    fname = dir_results + "scenario_" + hist.replace(" ", "_") + "_hi.pdf"
    plt.savefig(fname)
    plt.close("all")

    # utilities
    # ---

    plt.figure(figsize=(7, 3.8))
    ax = plt.subplot(111)
    ax.set_ylim((0, 0.8))  # NOTE hardcoded for ease of comparison

    # store each bar chart element to add legends later
    rects = list()

    # plot historical scenario
    UV = histD[hist]["UV"]
    colour = scenario2colour[hist]
    rect = ax.bar(x, UV, width=width, color=colour, align="center")
    rects.append(rect[0])

    # plot each voting-dynamical endpoint
    for idx, vot_sys in enumerate(voting_systemV):
        UV = UvotD[hist][vot_sys]
        colour = scenario2colour[vot_sys]
        rect = ax.bar(
            x + (idx + 1) * xdel, UV, width=width, color=colour, align="center"
        )
        rects.append(rect[0])

    # deocrate the plot and save
    ylabel = r"utilities $U_i$"
    xticks_w = x + xdel * (lv - 1) / 2
    xticklabels_w = [
        r"$w_" + str(foc + 1) + " = " + str(wV[foc]) + "$" for foc in range(n)
    ]
    ax.set_ylabel(ylabel, fontsize="x-large")
    ax.set_xticks(xticks_w)
    ax.set_xticklabels(xticklabels_w)
    ax.legend(rects, [hist] + voting_systemV)
    plt.tight_layout()
    fname = dir_results + "scenario_" + hist.replace(" ", "_") + "_ui.pdf"
    plt.savefig(fname)
    plt.close("all")
