# a figure to explain how the weightings affect willingness to trade off
# conservation versus economic outcomes

import matplotlib.pyplot as plt
import numpy as np

# ---
# parameters
# ---

# where the result directory is
dir_results = "../results/"

# conservation-interest parameters
w1 = 0.75
w2 = 0.25

# cost of production per unit effort
c = 0.05

# fixed utility
Ui = 0.5

# number of grid points for plotting
ngrid = 101


# ---
# plot explanatory figure
# ---

# x-axis "stock abundance" is (1 - H)
# y-axis "Net income" is h_i (1 - H - c)
# therefore Eq. 2 becomes
#   Ui = stock_abund^{w_i} * income^{1-w_i}

stock_abundV = np.linspace(0, 1, ngrid)
income_1V = [(Ui / (s**w1)) ** (1 / (1 - w1)) for s in stock_abundV]
income_2V = [(Ui / (s**w2)) ** (1 / (1 - w2)) for s in stock_abundV]

# plot it
plt.figure(figsize=(5, 3.8))
plt.plot(stock_abundV, income_1V, lw=2, color="red", ls="dashed", label=r"$w_1 = 0.75$")
plt.plot(stock_abundV, income_2V, lw=2, color="blue", label=r"$w_2 = 0.25$")
plt.ylim((0, 1.5))
plt.xlim((0, 1))
plt.legend(loc="best")
plt.xlabel(r"Stock abundance, $1-H$", fontsize="x-large")
plt.ylabel(r"Net income, $h_i (1 - H - c)$", fontsize="x-large")
plt.tight_layout()
plt.savefig(dir_results + "explain_weightings.pdf")
plt.close("all")
