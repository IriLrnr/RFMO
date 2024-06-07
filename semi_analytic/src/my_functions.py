# ---
# function for calculating the utilities
# ---


def UV_fnc(pars, hV):
    # unpack parameters
    wV = pars["wV"]
    c = pars["c"]

    # calculate the utility to each player given the hi vector
    H = sum(hV)
    UV = [(1 - H) ** wi * (hi * (1 - H - c)) ** (1 - wi) for wi, hi in zip(wV, hV)]

    return UV


# ---
# direct methods for solving the harvesting efforts (used for checking)
# ---


# Open-access derivatives, d U_i / d h_i, for Nash equilibrium
def dUdhiV_unconstrained_fnc(pars, hV):
    """
    The derivatives of each player's utility function with respect to their
    individual harvesting effort. When all of these are 0, we are at a singular
    strategy. This is used to find the Nash equilibrium.
    """

    # unpack parameters
    wV = pars["wV"]
    c = pars["c"]

    # total harvest effort is sum of harvesting efforts
    H = sum(hV)

    # each component function and its derivatives

    # conservation component
    AV = [(1 - H) ** wi for wi in wV]
    dAdhiV = [-wi * (1 - H) ** (wi - 1) for wi in wV]

    # final component of income
    BV = [(1 - c - H) ** (1 - wi) for wi in wV]
    dBdhiV = [-(1 - wi) * (1 - c - H) ** (-wi) for wi in wV]

    # middle component of income
    CV = [hi ** (1 - wi) for hi, wi in zip(hV, wV)]
    dCdhiV = [(1 - wi) * hi ** (-wi) for hi, wi in zip(hV, wV)]

    # product rule for derivatives
    dUdhiV = [
        dA * B * C + dB * A * C + dC * A * B
        for A, B, C, dA, dB, dC in zip(AV, BV, CV, dAdhiV, dBdhiV, dCdhiV)
    ]

    return dUdhiV


# Fixed-gamma derivative, d U_i / d H, for voting optima
def dUidH_fixed_gamma_fnc(pars, gammaV, i, H):
    """
    The derivative of focal player i's utility function with respect to the
    total harvesting effort where each player's effort is kept at fixed
    proportion of total effort gamma_i. If player i is the critical voter, then
    when i's derivative is zero, we have found a voting-dynamical singular
    total harvest rate \hat{H}.
    """

    # unpack parameters
    wi = pars["wV"][i]
    c = pars["c"]
    gammai = gammaV[i]

    # each component function and its derivatives

    # conservation component
    A = (1 - H) ** wi
    dAdH = -wi * (1 - H) ** (wi - 1)

    # final component of income
    B = (1 - c - H) ** (1 - wi)
    dBdH = -(1 - wi) * (1 - c - H) ** (-wi)

    # middle component of income
    C = (gammai * H) ** (1 - wi)
    dCdH = (1 - wi) * gammai ** (1 - wi) * H ** (-wi)

    # product rule
    dUidH = dAdH * B * C + dBdH * A * C + dCdH * A * B

    return dUidH


# ---
# more failsafe and direct methods for solving the harvesting efforts
# ---


def Hhat_fnc(pars, i):
    """
    Returns the left-root of the Q_HO_i function, which is player i's collective-harvesting optimum
    """

    # extract parameters
    c = pars["c"]
    wV = pars["wV"]

    left_root = (
        1
        / 2
        * (
            c
            + 2 * wV[i]
            + (4 * c * wV[i] ** 2 + c**2 - 8 * c * wV[i] + 2 * c + 1) ** (1 / 2)
            - 3
        )
        / (wV[i] - 2)
    )

    return left_root


def left_root_NE_fnc(pars, hV):
    """
    Find the zeros of this function to find the Nash equilibrium
    """

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
