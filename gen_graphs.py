"""Graphs for Light interception publication."""

import math

import matplotlib.pyplot as plt
import numpy as np
from environmental_biophysics.solar_radiation import (
    get_solar_radiation_interception_sub_daily,
)

from rad_competition_methods import (
    rad_intercpt_apsim,
    rad_intercpt_cycles,
    rad_intercpt_wallace,
)


def gen_fig1() -> None:
    """Figure 1."""
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
            math.e,
            3.11166667,
            3.5,
        ]
    )
    x_sp1 = 0.5
    x_sp2 = 2
    angles_deg = np.linspace(0, 90, 19)
    sp1, sp2 = get_solar_radiation_interception_sub_daily(
        atm_transm, atm_press, leaf_transm, leaf_area_index, x_sp1, x_sp2, angles_deg
    )
    tot_lai = 2 * leaf_area_index
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    apsim_sp1 = np.zeros(array_size)
    apsim_sp2 = np.zeros(array_size)
    for i in range(array_size):
        apsim_sp1[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[0]
        apsim_sp2[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[1]
    plt.plot(
        tot_lai,
        sp1,
        label=r"Sub-daily sp1 %s=0.5 $L$%%=%.0f" % (r"$\chi$", sp1_lai_percent * 100),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        tot_lai,
        sp2,
        label=r"Sub-daily sp2 %s=2 $L$%%=%.0f" % (r"$\chi$", sp1_lai_percent * 100),
        marker="v",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp1,
        label=r"Daily Cycles sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        apsim_sp2,
        label=r"Daily Cycles sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
    )
    plt.ylabel("Light interception", fontsize=14, labelpad=8)
    plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=14, labelpad=8)
    plt.xlim(0, 7)
    plt.ylim(0, 0.8)
    plt.legend(loc="upper left", prop={"size": 14}, frameon=False)
    plt.savefig("figures/Figure1.svg")


def gen_fig2() -> None:
    """Figure 2."""
    plt.figure(1, figsize=(8, 28))

    # Graph 2.1 : 50 / 50 lai; k1 = 0.4, k2 = 0.6
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    sp_height = np.ones(array_size)  # height equals 1
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    apsim_sp1 = np.zeros(array_size)
    apsim_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[0]
        apsim_sp1[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[0]
        wallace_sp2[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[1]
        apsim_sp2[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[1]
    plt.subplot(3, 1, 1).tick_params(axis="both", which="major", labelsize=20)
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp1,
        label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp2,
        label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
    )
    plt.ylabel("Light interception", fontsize=22, labelpad=8)
    plt.text(
        6, 0.9, "A", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22
    )
    plt.xlim(0, 7)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
    # Graph 2.2 : 80 / 20 lai; k1 = 0.4, k2 = 0.6
    sp1_lai_percent = 0.8
    sp2_lai_percent = 0.2
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    sp_height = np.ones(array_size)  # height equals 1
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    apsim_sp1 = np.zeros(array_size)
    apsim_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[0]
        apsim_sp1[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[0]
        wallace_sp2[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[1]
        apsim_sp2[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[1]
    plt.subplot(3, 1, 2).tick_params(axis="both", which="major", labelsize=20)
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp1,
        label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp2,
        label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
    )
    plt.ylabel("Light interception", fontsize=22, labelpad=8)
    plt.text(
        6, 0.8, "B", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22
    )
    plt.xlim(0, 7)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
    # Graph 1.3 : 20 / 80 lai; k1 = 0.4, k2 = 0.6
    sp1_lai_percent = 0.2
    sp2_lai_percent = 0.8
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    sp_height = np.ones(array_size)  # height equals 1
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    apsim_sp1 = np.zeros(array_size)
    apsim_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[0]
        apsim_sp1[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[0]
        wallace_sp2[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], sp_height[i]],
                [k_species2[i], lai_species2[i], sp_height[i]],
            ]
        )[1]
        apsim_sp2[i] = rad_intercpt_apsim(
            [[k_species1[i], lai_species1[i]], [k_species2[i], lai_species2[i]]]
        )[1]
    plt.subplot(3, 1, 3).tick_params(axis="both", which="major", labelsize=20)
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp1,
        label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f" % (k_sp1, sp1_lai_percent * 100),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        apsim_sp2,
        label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f" % (k_sp2, sp2_lai_percent * 100),
        marker="v",
        color="k",
    )
    plt.ylabel("Light interception", fontsize=22, labelpad=8)
    plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=22, labelpad=8)
    plt.text(
        6, 0.9, "C", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=22
    )
    plt.xlim(0, 7)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", prop={"size": 18}, frameon=False)
    plt.setp(plt.subplot(3, 1, 1).get_xticklabels(), visible=False)
    plt.setp(plt.subplot(3, 1, 2).get_xticklabels(), visible=False)
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("figures/Figure2.svg")


def gen_fig3() -> None:
    """Figure 3."""
    # Figure 3 Wallace with different heights
    plt.figure(2, figsize=(16, 12))
    # Graph 3.1: lai:50/50, k:0.5/0.5, sp_height:1/0.5
    plt.subplot(3, 3, 1).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.5
    k_sp2 = 0.5
    height_sp1 = 1
    height_sp2 = 0.5
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
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
    # Graph 3.2: lai:50/50, k:0.4/0.6, sp_height:0.5/1
    plt.subplot(3, 3, 2).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 0.5
    height_sp2 = 1
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.text(
        0.35, 0.7, "B", bbox={"facecolor": "white", "alpha": 0, "pad": 10}, fontsize=18
    )
    plt.xlim(0, 7)
    plt.ylim(0, 1.15)
    plt.legend(loc="upper left", prop={"size": 16}, frameon=False)

    # Graph 3.3: lai:50/50, k:0.4/0.6, sp_height:1/0.5
    plt.subplot(3, 3, 3).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 1
    height_sp2 = 0.5
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
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

    # Graph 3.4: lai:80/20, k:0.4/0.6, sp_height:0.5/1
    plt.subplot(3, 3, 4).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.8
    sp2_lai_percent = 0.2
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 0.5
    height_sp2 = 1
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
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

    # Graph 3.5: lai:20/80, k:0.4/0.6, sp_height:1/0.5
    plt.subplot(3, 3, 5).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.2
    sp2_lai_percent = 0.8
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 1
    height_sp2 = 0.5
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
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

    # Graph 3.6: lai:20/80, k:0.4/0.6, sp_height:0.5/1
    plt.subplot(3, 3, 6).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.2
    sp2_lai_percent = 0.8
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 0.5
    height_sp2 = 1
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"sp 1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"sp 2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
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


def gen_fig4() -> None:
    """Figure 4."""
    # Figure 4 Cycles and Wallace comparison
    plt.figure(3)
    plt.subplot(1, 1, 1).tick_params(axis="both", which="major", labelsize=16)
    sp1_lai_percent = 0.5
    sp2_lai_percent = 0.5
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2
    array_size = len(lai_species1)
    k_sp1 = 0.4
    k_sp2 = 0.6
    height_sp1 = 0.5
    height_sp2 = 1
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    wallace_sp1 = np.zeros(array_size)
    wallace_sp2 = np.zeros(array_size)
    cycles_sp1 = np.zeros(array_size)
    cycles_sp2 = np.zeros(array_size)
    for i in range(array_size):
        wallace_sp1[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        cycles_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[0]
        wallace_sp2[i] = rad_intercpt_wallace(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
        cycles_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
            ]
        )[1]
    plt.plot(
        lai_total,
        wallace_sp1,
        label=r"Wallace sp1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        cycles_sp1,
        label=r"Cycles sp1 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp1, sp1_lai_percent * 100, height_sp1),
        marker="o",
        color="k",
    )
    plt.plot(
        lai_total,
        wallace_sp2,
        label=r"Wallace sp2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
        marker="v",
        color="k",
        markerfacecolor="white",
    )
    plt.plot(
        lai_total,
        cycles_sp2,
        label=r"Cycles sp2 $k$=%.1f $L$%%=%.0f $h$=%.1f"
        % (k_sp2, sp2_lai_percent * 100, height_sp2),
        marker="v",
        color="k",
    )
    plt.ylabel("Light interception", fontsize=18, labelpad=8)
    plt.xlabel("Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=18, labelpad=8)
    plt.xlim(0, 7)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", prop={"size": 14}, frameon=False)
    plt.savefig("figures/Figure4.svg")


def gen_fig5() -> None:
    """Figure 5."""
    # Figure 5 Cycles and APSIM comparison
    fig4 = plt.figure(4)
    sp1_lai_percent = 0.33
    sp2_lai_percent = 0.33
    sp3_lai_percent = 0.33
    k_sp1 = 0.6
    k_sp2 = 0.6
    k_sp3 = 0.6
    height_sp1 = 1
    height_sp2 = 2
    height_sp3 = 4
    min_lai = 0.01
    max_lai = 7
    lai_species1 = np.linspace(
        min_lai * sp1_lai_percent, max_lai * sp1_lai_percent, num=10
    )
    lai_species2 = np.linspace(
        min_lai * sp2_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_species3 = np.linspace(
        min_lai * sp3_lai_percent, max_lai * sp2_lai_percent, num=10
    )
    lai_total = lai_species1 + lai_species2 + lai_species3
    array_size = len(lai_species1)
    k_species1 = np.ones(array_size) * k_sp1
    k_species2 = np.ones(array_size) * k_sp2
    k_species3 = np.ones(array_size) * k_sp3
    height_species1 = np.ones(array_size) * height_sp1
    height_species2 = np.ones(array_size) * height_sp2
    h_species3 = np.ones(array_size) * height_sp3
    apsim_sp1 = np.zeros(array_size)
    apsim_sp2 = np.zeros(array_size)
    apsim_sp3 = np.zeros(array_size)
    cycles_sp1 = np.zeros(array_size)
    cycles_sp2 = np.zeros(array_size)
    cycles_sp3 = np.zeros(array_size)
    for i in range(array_size):
        apsim_sp1[i] = rad_intercpt_apsim(
            [
                [k_species1[i], lai_species1[i]],
                [k_species2[i], lai_species2[i]],
                [k_species3[i], lai_species3[i]],
            ]
        )[0]
        apsim_sp2[i] = rad_intercpt_apsim(
            [
                [k_species1[i], lai_species1[i]],
                [k_species2[i], lai_species2[i]],
                [k_species3[i], lai_species3[i]],
            ]
        )[1]
        apsim_sp3[i] = rad_intercpt_apsim(
            [
                [k_species1[i], lai_species1[i]],
                [k_species2[i], lai_species2[i]],
                [k_species3[i], lai_species3[i]],
            ]
        )[2]
        cycles_sp1[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
                [k_species3[i], lai_species3[i], h_species3[i]],
            ]
        )[0]
        cycles_sp2[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
                [k_species3[i], lai_species3[i], h_species3[i]],
            ]
        )[1]
        cycles_sp3[i] = rad_intercpt_cycles(
            [
                [k_species1[i], lai_species1[i], height_species1[i]],
                [k_species2[i], lai_species2[i], height_species2[i]],
                [k_species3[i], lai_species3[i], h_species3[i]],
            ]
        )[2]
    axes1 = fig4.add_axes([0.1, 0.1, 0.8, 0.8])
    axes2 = fig4.add_axes([0.55, 0.55, 0.3, 0.3])
    axes1.plot(
        lai_total,
        apsim_sp1,
        label="APSIM (3 species)",
        marker="^",
        color="k",
        markerfacecolor="white",
    )
    axes1.plot(
        lai_total,
        cycles_sp1,
        label=r"Cycles sp 1 $h$=%.0f" % height_sp1,
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    axes1.plot(
        lai_total,
        cycles_sp2,
        label=r"Cycles sp 2 $h$=%.0f" % height_sp2,
        marker="o",
        color="k",
        markerfacecolor="gray",
    )
    axes1.plot(
        lai_total,
        cycles_sp3,
        label=r"Cycles sp 3 $h$=%.0f" % height_sp3,
        marker="o",
        color="k",
    )
    axes1.set_ylabel("Light interception", fontsize=16, labelpad=0)
    axes1.set_xlabel(
        "Total leaf area index " + r"(m$^2$ m$^{-2}$)", fontsize=16, labelpad=0
    )
    axes1.set_xlim([0, 7])
    axes1.set_ylim([0, 1])
    axes1.legend(loc="upper left", prop={"size": 12}, frameon=False)
    axes1.xaxis.set_tick_params(labelsize=12)
    axes1.yaxis.set_tick_params(labelsize=12)
    # Figure 4 Error between Cycles and APSIM comparison
    axes2.plot(
        lai_total,
        apsim_sp1 - apsim_sp1,
        label=r"APSIM (3 species)",
        marker="^",
        color="k",
        markerfacecolor="white",
    )
    axes2.plot(
        lai_total,
        cycles_sp1 - apsim_sp1,
        label=r"Cycles sp 1 $h$=%.0f" % height_sp1,
        marker="o",
        color="k",
        markerfacecolor="white",
    )
    axes2.plot(
        lai_total,
        cycles_sp2 - apsim_sp1,
        label=r"Cycles sp 2 $h$=%.0f" % height_sp2,
        marker="o",
        color="k",
        markerfacecolor="gray",
    )
    axes2.plot(
        lai_total,
        cycles_sp3 - apsim_sp1,
        label=r"Cycles sp 3 $h$=%.0f" % height_sp3,
        marker="o",
        color="k",
    )
    axes2.set_ylabel("Difference", fontsize=12, labelpad=8)
    axes2.set_yticks([-0.2, -0.1, 0.0, 0.1, 0.2])
    axes2.set_xticks([0, 1, 2, 3, 4, 5, 6, 7])
    axes2.yaxis.set_tick_params(labelsize=12)
    axes2.xaxis.set_tick_params(labelsize=12)
    plt.savefig("figures/Figure5.svg")


def get_fig6() -> None:
    """Figure 6."""
    # Barilot vs Cycles
    fig5 = plt.figure(5, figsize=(7, 9))
    axes1 = fig5.add_axes([0.1, 0.1, 0.8, 0.8])
    axes2 = fig5.add_axes([0.62, 0.13, 0.24, 0.20])
    fname = "data/barillot_data.csv"
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
    axes1.legend(loc="upper left", prop={"size": 14}, frameon=False)
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
    axes2.text(-30, 0.16, f"Mean pea bias = {pea_ave_abs_bias:.3f}", fontsize=11)
    axes2.text(-30, 0.18, f"Mean wheat bias = {wheat_ave_abs_bias:.3f}", fontsize=11)
    plt.savefig("figures/Figure6.svg")


if __name__ == "__main__":
    gen_fig1()
    gen_fig2()
    gen_fig3()
    gen_fig4()
    gen_fig5()
    get_fig6()
