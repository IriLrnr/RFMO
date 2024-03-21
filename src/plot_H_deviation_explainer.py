# use plots of the relationship between U_i(H) and H
# to explain how the voting dynamics works

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
# plot H curves and voting explainer
# ---

# create a grid with more points near the ends
# ---

H_lo = 0
H_hi = 1 - c
dgrid = 0.1
HV = (
    list(np.linspace(H_lo, H_lo + dgrid, ngrid))
    + list(np.linspace(H_lo + dgrid, H_hi - dgrid, ngrid))
    + list(np.linspace(H_hi - dgrid, H_hi, ngrid))
)

# controls for plots
# ---

prop_cycle = plt.rcParams['axes.prop_cycle']
colourV = prop_cycle.by_key()['color']
# xlim = (-0.05, H_hi+0.05)
xlim = (0, H_hi)
ylim_util = (0, 0.9)
ylim_vote_count = (-0.3, 7.8)

# count how many votes from left and right
vote_count_HV = [1-c]
vote_count_leftV = [n]
vote_count_rightV = [0]
for i, Hhat in enumerate(HhatV):
    vote_count_HV += [Hhat, Hhat]
    vote_count_leftV += [n-i, n-i-1]
    vote_count_rightV += [i, i+1]
vote_count_HV += [0]
vote_count_leftV += [0]
vote_count_rightV += [n]


# stack the data
gammaM = [gamma_EXV, gamma_OAV]
titleV = ["historically exploratory", "historically open access"]
vote_countM = [vote_count_rightV, vote_count_leftV]
H_hist_labelV = [r"$H^*$", r"$H^{(0)}$"]

fig=plt.figure()
ax_tl = plt.subplot(221)
ax_tr = plt.subplot(222, sharex = ax_tl)
ax_bl = plt.subplot(223)
ax_br = plt.subplot(224)
axs = [[ax_tl, ax_tr], [ax_bl, ax_br]]

# annotate historical and vote effect
# ---
H0 = n * 1e-3

axs[0][0].axvline(H0, ls="dotted", color="gray")
#axs[0][0].text(H0 - 0.01, ylim_util[1], r"historical, $H^{(0)}$ ", ha="right", va="top", fontsize="xx-small", rotation="vertical")
axs[0][0].text(H0, 0.8, r"$\longrightarrow$", ha="left", va="center", fontsize="x-large")
axs[0][0].text(0.15, 0.8, "voting increases harvest", ha="left", va="center", fontsize="small")

axs[0][1].axvline(Hs, ls="dotted", color="gray")
# axs[0][1].text(Hs + 0.01, ylim_util[1], r"historical, $H^{*}$ ", ha="left", va="top", fontsize="xx-small", rotation="vertical")
axs[0][1].text(Hs - 0.01, 0.7, r"historical ", ha="right", va="top", fontsize="small", rotation="vertical")
axs[0][1].text(Hs, 0.8, r"$\longleftarrow$", ha="right", va="center", fontsize="x-large")
axs[0][1].text(Hs - 0.15, 0.8, "voting decreases harvest", ha="right", va="center", fontsize="small")

# axis labels
ax_tl.set_ylabel("utility to each country")
ax_tl.set_yticks([0, 0.9])
ax_tr.set_yticks([0, 1], ["", ""])
ax_tl.set_xticks([0, H_hi], ['', ''])

ax_bl.set_yticks([0, 4, 6, 7])
ax_br.set_yticks([0, 4, 6, 7], [""]*4)
ax_bl.set_ylabel("nbr voting for change")

# plot curves and vote count
# ---

for scenario in [0, 1]:
    
    # extract needed items for this scenario
    gammaV = gammaM[scenario]
    vote_countV = vote_countM[scenario]
    title = titleV[scenario]


    # get the curve of utilities versus H
    # ---

    UM = [None] * ngrid * 3
    for idx, H_new in enumerate(HV):
        # the new hi's will be h_i = gamma_i * H_new
        hV = [gamma * H_new for gamma in gammaV]

        # vector of utilities (u_1, ..., u_n)
        UV = UV_fnc(pars, hV)

        # store
        UM[idx] = UV

    # plot each country's utility curve
    # ---

    axs[0][scenario].set_title(title)
    UNM = np.array(UM)
    for foc in range(n):
        # how focal player's utility varies with H
        UV = UNM[:, foc]
        colour = colourV[foc]

        # plot how utility varies with collective-harvesting effort
        axs[0][scenario].plot(HV, UV) #, label=r"$i = " + str(foc + 1) + r"$", color=colour)

        # plot a dot at focal player's optimum
        Hhat = HhatV[foc]
        wi = wV[foc]
        hi = gammaV[foc] * Hhat
        U_Hhat = (1 - Hhat) ** wi * (hi * (1 - Hhat - c)) ** (1 - wi)
        axs[0][scenario].scatter([Hhat], [U_Hhat], s=10, color=colour)
        axs[0][scenario].plot([Hhat, Hhat], [0, U_Hhat], alpha=0.2, color=colour)


    # plot the vote counts, mark where consensus etc is
    # ---

    # stepping line showing how many vote for change
    # axs[1][scenario].plot(HV, vote_countV, color="black")
    axs[1][scenario].set_xlabel("harvest effort $H$")
    axs[1][scenario].plot(vote_count_HV, vote_countV, color="black")

    # light lines showing H hat 
    for foc in range(n):
        # axs[1][scenario].axvline(HhatV[foc], alpha=0.2, color=colourV[foc])
        axs[1][0].plot([HhatV[foc], HhatV[foc]], [foc, 10], alpha=0.2, color=colourV[foc])
        axs[1][1].plot([HhatV[foc], HhatV[foc]], [n-foc, 10], alpha=0.2, color=colourV[foc])


    # lines at consensus, supermajority, majority
    # ---
    
    # to be used for indicating voting 
    colour = "black"
    ls = "dotted"

    # left panel above
    axs[1][0].plot([HhatV[-1], 10], [7, 7], ls=ls, color=colour)
    axs[1][0].plot([HhatV[-2], 10], [6, 6], ls=ls, color=colour)
    axs[1][0].plot([HhatV[-4], 10], [4, 4], ls=ls, color=colour)

    # right panel above
    axs[1][1].plot([-10, HhatV[0]], [7, 7], ls=ls, color=colour)
    axs[1][1].plot([-10, HhatV[1]], [6, 6], ls=ls, color=colour)
    axs[1][1].plot([-10, HhatV[3]], [4, 4], ls=ls, color=colour)

    # indicate where consensus, supermjority, are in left panel
    axs[1][0].text(H_hi-0.07, 7, "consensus", backgroundcolor="white", va="center", ha="right", fontsize="small")
    axs[1][0].text(H_hi-0.07, 6, "supermajority", backgroundcolor="white", va="center", ha="right", fontsize="small")
    axs[1][0].text(H_hi-0.07, 4, "majority", backgroundcolor="white", va="center", ha="right", fontsize="small")

    # voting-dynamical endpoints extending down - left panel
    axs[1][0].plot([HhatV[-1], HhatV[-1]], [-10, 7], ls=ls, color=colourV[n-1])
    axs[1][0].plot([HhatV[-2], HhatV[-2]], [-10, 6], ls=ls, color=colourV[n-2])
    axs[1][0].plot([HhatV[-4], HhatV[-4]], [-10, 4], ls=ls, color=colourV[n-4])
    ax_bl.set_xticks(
        [0, H0, HhatV[-1], HhatV[-2], HhatV[-4], H_hi], 
        labels=["", r"$H^{(0)}$ ", r"     $H_C^{\bullet}$", r"   $H_S^{\bullet}$", r"    $H_M^{\bullet}$", "1-c"]
    )

    # voting-dynamical endpoints extending down - right panel
    axs[1][1].plot([HhatV[0], HhatV[0]], [-10, 7], ls=ls, color=colourV[0])
    axs[1][1].plot([HhatV[1], HhatV[1]], [-10, 6], ls=ls, color=colourV[1])
    axs[1][1].plot([HhatV[3], HhatV[3]], [-10, 4], ls=ls, color=colourV[3])
    ax_br.set_xticks(
        [0, HhatV[0], HhatV[1], HhatV[3], Hs, H_hi], 
        labels=["0", r"      $H_C^{\bullet}$", r"$H_S^{\bullet}$", r"$H_M^{\bullet}$", r"$H^*$", "1-c"]
    )

    # put in arrows
    axs[1][0].text(H0, 7.3, r"$\longrightarrow$", ha="left", va="center", fontsize="x-large")
    axs[1][1].text(Hs, 7.3, r"$\longleftarrow$", ha="right", va="center", fontsize="x-large")

    # historical harvest coming down
    axs[1][1].plot([Hs, Hs], [7, 10], ls=ls, color="gray")
    axs[1][0].plot([H0, H0,], [7, 10], ls=ls, color="gray")

    # axis limits
    # ---

    axs[0][scenario].set_xlim(xlim)
    axs[0][scenario].set_ylim(ylim_util)
    axs[1][scenario].set_xlim(xlim)
    axs[1][scenario].set_ylim(ylim_vote_count)


# partial legend in the first panel
# ---

axs[0][0].plot([], [], color=colourV[0], label=r"$U_1$")
axs[0][0].plot([], [], color=colourV[1], label=r"$U_2$")
axs[0][0].plot([], [], color="white", label=r"$\vdots$")
axs[0][0].plot([], [], color=colourV[6], label=r"$U_7$")
axs[0][0].scatter([], [], s=10, color=colourV[0], label=r"$\hat{H}_1$")
axs[0][0].plot([], [], color="white", label=r"$\vdots$")
axs[0][0].scatter([], [], s=10, color=colourV[5], label=r"$\hat{H}_6$")
axs[0][0].scatter([], [], s=10, color=colourV[6], label=r"$\hat{H}_7$")
axs[0][0].legend(loc="best", fontsize="x-small")



plt.tight_layout()
plt.savefig(dir_results + "H_deviation_explainer.pdf")
plt.close("all")
