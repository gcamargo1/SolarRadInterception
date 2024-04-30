"""Graphs for Light interception publication"""
import matplotlib.pyplot as plt
import numpy as np

from rad_competition_methods import (
    rad_intercpt_apsim,
    rad_intercpt_wallace,
    rad_intercpt_cycles,
)
from environ_biophysics import rad_intercpt_sub_daily

# Figure 0
plt.figure(0, figsize=[8, 8])
plt.subplot(1, 1, 1).tick_params(axis="both", which="major", labelsize=12)

atm_transm = 0.75
atm_press = 101.3
leaf_transm = 0.8

leaf_area_index = np.array(
    [
        0.005,
        0.39333333,
        0.78166667,
        1.17,
        1.55833333,
        1.94666667,
        2.335,
        2.72333333,
        3.11166667,
        3.5,
    ]
)
X_SP1 = 0.5
X_SP2 = 2
angles_deg = np.linspace(0, 90, 19)

SP1, SP2 = rad_intercpt_sub_daily(
    atm_transm, atm_press, leaf_transm, leaf_area_index, X_SP1, X_SP2, angles_deg
)

TOT_LAI = 2 * leaf_area_index
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
SP_HEIGHT = np.ones(ARRAY_SIZE)  # height equals 1
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
APSIM_SP1 = np.zeros(ARRAY_SIZE)
APSIM_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    APSIM_SP1[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[0]
    APSIM_SP2[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[1]
plt.plot(
    TOT_LAI,
    SP1,
    label=r"Sub-daily sp1 %s=0.5 $L$%%=%.0f" % (r"$\chi$", SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    TOT_LAI,
    SP2,
    label=r"Sub-daily sp2 %s=2 $L$%%=%.0f" % (r"$\chi$", SP1_LAI_PERCENT * 100),
    marker="v",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP1,
    label=r"Daily Cycles sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP2,
    label=r"Daily Cycles sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
)
plt.ylabel("Light interception", fontsize=14, labelpad=8)
plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=14, labelpad=8)
plt.xlim(0, 7)
plt.ylim(0, 0.8)
plt.legend(loc="upper left", prop={"size": 14}, frameon=False)
plt.savefig("figures/Figure1.svg")
# Figure 1
plt.figure(1, figsize=(8, 28))
# Graph 1.1 : 50 / 50 lai; k1 = 0.4, k2 = 0.6
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
SP_HEIGHT = np.ones(ARRAY_SIZE)  # height equals 1
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
APSIM_SP1 = np.zeros(ARRAY_SIZE)
APSIM_SP2 = np.zeros(ARRAY_SIZE)

for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[0]
    APSIM_SP1[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[1]
    APSIM_SP2[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[1]

plt.subplot(3, 1, 1).tick_params(axis="both", which="major", labelsize=20)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP1,
    label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP2,
    label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
)
plt.ylabel("Light interception", fontsize=22, labelpad=8)
# plt.xlabel('Total leaf area index ' + r'(m$^2$ m$^{-2}$)', fontsize=15,
#           labelpad=8)
plt.text(6, 0.9, "A", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22)
plt.xlim(0, 7)
plt.ylim(0, 1)
plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
# Graph 1.2 : 80 / 20 lai; k1 = 0.4, k2 = 0.6
SP1_LAI_PERCENT = 0.8
SP2_LAI_PERCENT = 0.2
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
SP_HEIGHT = np.ones(ARRAY_SIZE)  # height equals 1
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
APSIM_SP1 = np.zeros(ARRAY_SIZE)
APSIM_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[0]
    APSIM_SP1[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[1]
    APSIM_SP2[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[1]


plt.subplot(3, 1, 2).tick_params(axis="both", which="major", labelsize=20)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP1,
    label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP2,
    label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
)
plt.ylabel("Light interception", fontsize=22, labelpad=8)
# plt.xlabel('Total leaf area index ' + r'(m$^2$ m$^{-2}$)', fontsize=15,
#           labelpad=8)
plt.text(6, 0.8, "B", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22)
plt.xlim(0, 7)
plt.ylim(0, 1)
plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
# Graph 1.3 : 20 / 80 lai; k1 = 0.4, k2 = 0.6
SP1_LAI_PERCENT = 0.2
SP2_LAI_PERCENT = 0.8
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
SP_HEIGHT = np.ones(ARRAY_SIZE)  # height equals 1
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
APSIM_SP1 = np.zeros(ARRAY_SIZE)
APSIM_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[0]
    APSIM_SP1[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], SP_HEIGHT[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], SP_HEIGHT[i]],
        ]
    )[1]
    APSIM_SP2[i] = rad_intercpt_apsim(
        [[K_SPECIES1[i], LAI_SPECIES1[i]], [K_SPECIES2[i], LAI_SPECIES2[i]]]
    )[1]
plt.subplot(3, 1, 3).tick_params(axis="both", which="major", labelsize=20)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP1,
    label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (K_SP1, SP1_LAI_PERCENT * 100),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    APSIM_SP2,
    label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (K_SP2, SP2_LAI_PERCENT * 100),
    marker="v",
    color="k",
)
plt.ylabel("Light interception", fontsize=22, labelpad=8)
plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=22, labelpad=8)
plt.text(6, 0.9, "C", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22)
plt.xlim(0, 7)
plt.ylim(0, 1)
plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
plt.setp(plt.subplot(3, 1, 1).get_xticklabels(), visible=False)
plt.setp(plt.subplot(3, 1, 2).get_xticklabels(), visible=False)
plt.subplots_adjust(hspace=0.05)
plt.savefig("figures/Figure2.svg")


# Figure 2 Wallace with different heights
plt.figure(2, figsize=(16, 12))
# Graph 2.1: lai:50/50, k:0.5/0.5, SP_HEIGHT:1/0.5
plt.subplot(3, 3, 1).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.5
K_SP2 = 0.5
HEIGHT_SP1 = 1
HEIGHT_SP2 = 0.5
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.text(
    0.35, 0.7, "A", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
plt.ylabel("Light interception", fontsize=18, labelpad=8)
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Graph 2.2: lai:50/50, k:0.4/0.6, SP_HEIGHT:0.5/1
plt.subplot(3, 3, 2).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 0.5
HEIGHT_SP2 = 1
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.text(
    0.35, 0.7, "B", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
# plt.ylabel('Light interception')
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Graph 2.3: lai:50/50, k:0.4/0.6, SP_HEIGHT:1/0.5
plt.subplot(3, 3, 3).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 1
HEIGHT_SP2 = 0.5
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.text(
    0.35, 0.7, "C", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Graph 2.4: lai:80/20, k:0.4/0.6, SP_HEIGHT:0.5/1
plt.subplot(3, 3, 4).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.8
SP2_LAI_PERCENT = 0.2
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 0.5
HEIGHT_SP2 = 1
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8)
plt.ylabel("Light interception", fontsize=18, labelpad=8)
plt.text(
    0.35, 0.7, "D", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Graph 2.5: lai:20/80, k:0.4/0.6, SP_HEIGHT:1/0.5
plt.subplot(3, 3, 5).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.2
SP2_LAI_PERCENT = 0.8
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 1
HEIGHT_SP2 = 0.5
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)

for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)

plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8)
plt.text(
    0.35, 0.7, "E", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Graph 2.6: lai:20/80, k:0.4/0.6, SP_HEIGHT:0.5/1
plt.subplot(3, 3, 6).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.2
SP2_LAI_PERCENT = 0.8
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 0.5
HEIGHT_SP2 = 1
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)

for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="white",
)

plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8)
plt.text(
    0.35, 0.7, "F", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
)
plt.xlim(0, 7)
plt.ylim(0, 1.15)
plt.legend(loc="upper left", prop={"size": 16}, frameon=False)
# Axis removal controls
plt.setp(plt.subplot(3, 3, 1).get_xticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 2).get_yticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 3).get_yticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 3).get_xticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 5).get_yticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 2).get_xticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 2).get_yticklabels(), visible=False)
plt.setp(plt.subplot(3, 3, 6).get_yticklabels(), visible=False)
plt.subplots_adjust(wspace=0.07, hspace=0.1)
plt.savefig("figures/Figure3.svg")
# Figure 3 Cycles and Wallace comparison
plt.figure(3)
plt.subplot(1, 1, 1).tick_params(axis="both", which="major", labelsize=16)
SP1_LAI_PERCENT = 0.5
SP2_LAI_PERCENT = 0.5
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2
ARRAY_SIZE = len(LAI_SPECIES1)
K_SP1 = 0.4
K_SP2 = 0.6
HEIGHT_SP1 = 0.5
HEIGHT_SP2 = 1
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
WALLACE_SP1 = np.zeros(ARRAY_SIZE)
WALLACE_SP2 = np.zeros(ARRAY_SIZE)
CYCLES_SP1 = np.zeros(ARRAY_SIZE)
CYCLES_SP2 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    WALLACE_SP1[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    CYCLES_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[0]
    WALLACE_SP2[i] = rad_intercpt_wallace(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
    CYCLES_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
        ]
    )[1]
plt.plot(
    LAI_TOTAL,
    WALLACE_SP1,
    label=r"Wallace SP1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    CYCLES_SP1,
    label=r"Cycles SP1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP1, SP1_LAI_PERCENT * 100, HEIGHT_SP1),
    marker="o",
    color="k",
)
plt.plot(
    LAI_TOTAL,
    WALLACE_SP2,
    label=r"Wallace SP2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="v",
    color="k",
    markerfacecolor="white",
)
plt.plot(
    LAI_TOTAL,
    CYCLES_SP2,
    label=r"Cycles SP2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
    % (K_SP2, SP2_LAI_PERCENT * 100, HEIGHT_SP2),
    marker="v",
    color="k",
)
plt.ylabel("Light interception", fontsize=18, labelpad=8)
plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8)
plt.xlim(0, 7)
plt.ylim(0, 1)
plt.legend(loc="upper left", prop={"size": 14}, frameon=False)
plt.savefig("figures/Figure4.svg")
# Figure 4 Cycles and APSIM comparison
fig4 = plt.figure(4)
SP1_LAI_PERCENT = 0.33
SP2_LAI_PERCENT = 0.33
SP3_LAI_PERCENT = 0.33
K_SP1 = 0.6
K_SP2 = 0.6
K_SP3 = 0.6
HEIGHT_SP1 = 1
HEIGHT_SP2 = 2
HEIGHT_SP3 = 4
MIN_LAI = 0.01
MAX_LAI = 7
LAI_SPECIES1 = np.linspace(MIN_LAI * SP1_LAI_PERCENT, MAX_LAI * SP1_LAI_PERCENT, num=10)
LAI_SPECIES2 = np.linspace(MIN_LAI * SP2_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_SPECIES3 = np.linspace(MIN_LAI * SP3_LAI_PERCENT, MAX_LAI * SP2_LAI_PERCENT, num=10)
LAI_TOTAL = LAI_SPECIES1 + LAI_SPECIES2 + LAI_SPECIES3
ARRAY_SIZE = len(LAI_SPECIES1)
K_SPECIES1 = np.ones(ARRAY_SIZE) * K_SP1
K_SPECIES2 = np.ones(ARRAY_SIZE) * K_SP2
K_SPECIES3 = np.ones(ARRAY_SIZE) * K_SP3
HEIGHT_SPECIES1 = np.ones(ARRAY_SIZE) * HEIGHT_SP1
HEIGHT_SPECIES2 = np.ones(ARRAY_SIZE) * HEIGHT_SP2
H_SPECIES3 = np.ones(ARRAY_SIZE) * HEIGHT_SP3
APSIM_SP1 = np.zeros(ARRAY_SIZE)
APSIM_SP2 = np.zeros(ARRAY_SIZE)
APSIM_SP3 = np.zeros(ARRAY_SIZE)
CYCLES_SP1 = np.zeros(ARRAY_SIZE)
CYCLES_SP2 = np.zeros(ARRAY_SIZE)
CYCLES_SP3 = np.zeros(ARRAY_SIZE)
for i in range(ARRAY_SIZE):
    APSIM_SP1[i] = rad_intercpt_apsim(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i]],
        ]
    )[0]
    APSIM_SP2[i] = rad_intercpt_apsim(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i]],
        ]
    )[1]
    APSIM_SP3[i] = rad_intercpt_apsim(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i]],
        ]
    )[2]
    CYCLES_SP1[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i], H_SPECIES3[i]],
        ]
    )[0]
    CYCLES_SP2[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i], H_SPECIES3[i]],
        ]
    )[1]
    CYCLES_SP3[i] = rad_intercpt_cycles(
        [
            [K_SPECIES1[i], LAI_SPECIES1[i], HEIGHT_SPECIES1[i]],
            [K_SPECIES2[i], LAI_SPECIES2[i], HEIGHT_SPECIES2[i]],
            [K_SPECIES3[i], LAI_SPECIES3[i], H_SPECIES3[i]],
        ]
    )[2]

axes1 = fig4.add_axes([0.1, 0.1, 0.8, 0.8])
axes2 = fig4.add_axes([0.55, 0.55, 0.3, 0.3])

axes1.plot(
    LAI_TOTAL,
    APSIM_SP1,
    label="APSIM (3 species)",
    marker="^",
    color="k",
    markerfacecolor="white",
)
axes1.plot(
    LAI_TOTAL,
    CYCLES_SP1,
    label=r"Cycles sp 1 $h$=%.0f" % (HEIGHT_SP1),
    marker="o",
    color="k",
    markerfacecolor="white",
)
axes1.plot(
    LAI_TOTAL,
    CYCLES_SP2,
    label=r"Cycles sp 2 $h$=%.0f" % (HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="gray",
)
axes1.plot(
    LAI_TOTAL,
    CYCLES_SP3,
    label=r"Cycles sp 3 $h$=%.0f" % (HEIGHT_SP3),
    marker="o",
    color="k",
)
axes1.set_ylabel("Light interception", fontsize=18, labelpad=8)
axes1.set_xlabel(
    "Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8
)
axes1.set_xlim([0, 7])
axes1.set_ylim([0, 1])
axes1.legend(loc="upper left", prop={"size": 16}, frameon=False)
axes1.xaxis.set_tick_params(labelsize=16)
axes1.yaxis.set_tick_params(labelsize=16)
# Figure 4 Error between Cycles and APSIM comparison
axes2.plot(
    LAI_TOTAL,
    APSIM_SP1 - APSIM_SP1,
    label=r"APSIM (3 species)",
    marker="^",
    color="k",
    markerfacecolor="white",
)
axes2.plot(
    LAI_TOTAL,
    CYCLES_SP1 - APSIM_SP1,
    label=r"Cycles sp 1 $h$=%.0f" % (HEIGHT_SP1),
    marker="o",
    color="k",
    markerfacecolor="white",
)
axes2.plot(
    LAI_TOTAL,
    CYCLES_SP2 - APSIM_SP1,
    label=r"Cycles sp 2 $h$=%.0f" % (HEIGHT_SP2),
    marker="o",
    color="k",
    markerfacecolor="gray",
)
axes2.plot(
    LAI_TOTAL,
    CYCLES_SP3 - APSIM_SP1,
    label=r"Cycles sp 3 $h$=%.0f" % (HEIGHT_SP3),
    marker="o",
    color="k",
)
axes2.set_ylabel("Difference", fontsize=18, labelpad=8)
axes2.set_yticks([-0.2, -0.1, 0.0, 0.1, 0.2])
axes2.yaxis.set_tick_params(labelsize=16)
axes2.xaxis.set_tick_params(labelsize=16)
plt.savefig("figures/Figure5.svg")


# Barilot vs Cycles
fig5 = plt.figure(5, figsize=(7, 9))
axes1 = fig5.add_axes([0.1, 0.1, 0.8, 0.8])
axes2 = fig5.add_axes([0.62, 0.13, 0.24, 0.20])
fname = "barillot_data.csv"

data = np.loadtxt(fname, dtype="float", delimiter=",", skiprows=1)

therm_time = data[:, 0]
pea_height = data[:, 1]
wheat_height = data[:, 2]
pea_light_intercpt = data[:, 3]
wheat_light_intercpt = data[:, 4]
pea_lai = data[:, 5]
wheat_lai = data[:, 6]

pea_sim_li = np.zeros(len(therm_time))
wheat_sim_li = np.zeros(len(therm_time))
length = range(len(therm_time))
kp_final = 0.540909090909  # values obtained through optimization
kw_final = 0.510606060606
for i in length:
    pea_sim_li[i], wheat_sim_li[i] = rad_intercpt_cycles(
        (
            [kp_final, pea_lai[i], pea_height[i]],
            [kw_final, wheat_lai[i], wheat_height[i]],
        )
    )

axes1.plot(
    therm_time,
    pea_sim_li + wheat_sim_li,
    label="Mixture Cycles",
    color="k",
    dashes=(5, 5),
    linewidth=2,
)
axes1.plot(
    therm_time,
    pea_light_intercpt + wheat_light_intercpt,
    label="Mixture Barillot",
    color="k",
    marker="s",
    markerfacecolor="white",
    markersize=5,
    linewidth=0,
)
axes1.plot(
    therm_time,
    pea_sim_li,
    label="Pea Cycles",
    color="k",
    dashes=[5, 3, 1, 3],
    linewidth=2,
)
axes1.plot(
    therm_time,
    pea_light_intercpt,
    label="Pea Barillot",
    color="k",
    marker="o",
    markerfacecolor="white",
    markersize=4,
    linewidth=0,
)
axes1.plot(
    therm_time,
    wheat_sim_li,
    label="Wheat Cycles",
    color="k",
    dashes=[1, 3],
    linewidth=2,
)
axes1.plot(
    therm_time,
    wheat_light_intercpt,
    label="Wheat Barillot",
    color="k",
    marker="^",
    markerfacecolor="k",
    markersize=4,
    linewidth=0,
)

axes1.set_xlim(0, 2000)
axes1.set_ylim(0, 1)
axes1.xaxis.set_tick_params(labelsize=16)
axes1.yaxis.set_tick_params(labelsize=16)
axes1.set_xlabel("Thermal time (C-day)", fontsize=18)
axes1.set_ylabel("Light interception", fontsize=18)
# axes1.subplots_adjust(hspace=0.1)
axes1.legend(loc="upper left", prop={"size": 14}, frameon=False)

total_diff = (pea_sim_li + wheat_sim_li) - (pea_light_intercpt + wheat_light_intercpt)
pea_diff = abs(pea_sim_li - pea_light_intercpt)
wheat_diff = abs(wheat_sim_li - wheat_light_intercpt)
pea_ave_abs_bias = pea_diff.mean()
wheat_ave_abs_bias = wheat_diff.mean()

axes2.plot(
    therm_time,
    pea_diff,
    label="Pea",
    color="k",
    marker="8",
    markerfacecolor="white",
    markersize=3,
    linewidth=0,
)
axes2.plot(
    therm_time,
    wheat_diff,
    ":",
    label="Wheat",
    color="k",
    marker="^",
    markersize=3,
    linewidth=0.1,
)
axes2.legend(loc="upper left", prop={"size": 11}, frameon=False)
axes2.set_yticks([0.0, 0.05, 0.1, 0.15])
axes2.set_xticks([0, 1000, 2000])
axes2.set_ylabel("Absolute difference", fontsize=15, labelpad=8)
axes2.text(-30, 0.16, "Mean pea bias = %.3f" % pea_ave_abs_bias, fontsize=11)
axes2.text(-30, 0.18, "Mean wheat bias = %.3f" % wheat_ave_abs_bias, fontsize=11)

plt.savefig("figures/Figure6.svg")

# plt.show()
