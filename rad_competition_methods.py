"""Solar radiation interception methods."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from environmental_biophysics.solar_radiation import get_height_weight_factor

if TYPE_CHECKING:
    from numpy.typing import NDArray


def rad_intercpt_cycles(
    crop_list: tuple[
        list[float | int | NDArray], list[float | int | NDArray], list[float | NDArray]
    ],
) -> NDArray:
    """Returns solar radiation intercepted on each species.

    Args:
        crop_list: list of species characteristics with three items:
            radiation extinction coefficients
            leaf area index [m2/m2]
            plant height [m]

    References:
         Camargo, G.G.T. 2015. Ph.D. Dissertation. Penn State University.
          https://etda.libraries.psu.edu/files/final_submissions/10226

    Examples:
        >>> rad_intercpt_cycles(([0.5,1,1],[0.5,1,2],[0.5,1,1]))
        array([0.23758411, 0.30170163, 0.23758411])
        >>> rad_intercpt_cycles(([0.5,1,0.5],[0.6,1.2,1],[0.7,1.4,1.5]))
        array([0.13822214, 0.29363746, 0.45733724])
    """
    # Variables init
    number_species = len(crop_list)
    extinction_coeff = np.zeros(number_species)
    leaf_area_index = np.zeros(number_species)
    height = np.zeros(number_species)
    transm_rad = np.zeros(number_species)  # if plant was alone
    rad_intercpt_dom = np.zeros(number_species)  # dominant species
    rad_intercpt_suppr = np.zeros(number_species)  # supressed species
    rad_intercpt = np.zeros(number_species)  # species rad interception
    transm_rad_others = np.ones(number_species)  # non-domnt sp rad transm
    height_dom = np.zeros(number_species)  # species canopy dominance factor
    dominant_fact = np.zeros(number_species)  # Dominant weight factor
    suppressed_fact = np.zeros(number_species)  # Suppressed weight factor
    hght_wght_fct = np.zeros(number_species)
    hght_wght_fct_adj = np.zeros(number_species)
    k_lai_prod = np.zeros(number_species)
    k_lai_prod_adj = np.zeros(number_species)
    total_transm = 1

    # Read and store inputs
    for i in range(number_species):
        extinction_coeff[i] = crop_list[i][0]
        leaf_area_index[i] = crop_list[i][1]
        height[i] = crop_list[i][2]

        # Transmitted radiation if all species had same height
        transm_rad[i] = math.exp(-extinction_coeff[i] * leaf_area_index[i])

        # Intercepted radiation if species was dominant
        rad_intercpt_dom[i] = 1 - transm_rad[i]
        k_lai_prod[i] = extinction_coeff[i] * leaf_area_index[i]

    # Calculate total transmitance, interception and height dominance
    for i in range(number_species):
        # Total transmitance if all species had the same height
        total_transm *= transm_rad[i]

        # Height dominance factor
        height_dom[i] = number_species * height[i] / height.sum()

    # Total radiation interception if species had the same height
    total_interception = 1 - total_transm

    # All species but ith species rad transmission
    for i in range(number_species):
        all_species = list(range(number_species))  # list with all species
        all_species.pop(i)  # remove ith species from list
        all_species_but_ith = all_species  # list of non-dominant species
        total_transm_but_ith_sp = 1  # sum of non-dominant species transmission
        for j in all_species_but_ith:
            total_transm_but_ith_sp *= transm_rad[j]

        # Total transmitted radiation from all species but ith
        transm_rad_others[i] = total_transm_but_ith_sp

    # Radiation interception by suppressed species once all other species
    # intercepts the radiation first
    for i in range(number_species):
        rad_intercpt_suppr[i] = rad_intercpt_dom[i] * transm_rad_others[i]

    # Determine two extremes weighing factors: in dominant_fact species will
    # intercept all radiation than it can based on k and LAI, in
    # suppressed_factor, species will only intercept radiation after all the
    # other species intercepted all rad that was possible
    for i in range(number_species):
        dominant_fact[i] = (
            rad_intercpt_dom[i] / total_interception * k_lai_prod.sum() / k_lai_prod[i]
        )
        suppressed_fact[i] = (
            rad_intercpt_suppr[i]
            / total_interception
            * k_lai_prod.sum()
            / k_lai_prod[i]
        )

        # Based on species height determine a height weight factor in between
        # dominant_fact and suppressed_fact values usin linear interpolation
        hght_wght_fct[i] = get_height_weight_factor(
            height_dom[i], dominant_fact[i], suppressed_fact[i], number_species
        )
        # Adjust extinction coefficient and leaf area index product
        k_lai_prod_adj[i] = extinction_coeff[i] * leaf_area_index[i] * hght_wght_fct[i]

    for i in range(number_species):
        # Adjust height weighting factor
        hght_wght_fct_adj[i] = (
            hght_wght_fct[i] / k_lai_prod_adj.sum() * k_lai_prod.sum()
        )

        # Radiation interception for each species
        rad_intercpt[i] = (
            total_interception * hght_wght_fct_adj[i] * k_lai_prod[i] / k_lai_prod.sum()
        )
    return rad_intercpt


def rad_intercpt_wallace(
    crop_list: tuple[list[float | int], list[float | int]],
) -> list[float]:
    """Returns radiation intercepted on each species.

    Args:
        crop_list: list of species characteristics:
            extinction_coeff: rad extinction coefficient
            leaf_area_index: leaf area index [m2/m2]
            height: plant height [m]

    References:
        Wallace, J.S., 1997. Evaporation and radiation interception by neighbouring
         plants. Quarterly Journal of the Royal Meteorological Society 123, 1885-1905.

    Examples:
        >>> rad_intercpt_wallace(([0.5, 1, 1], [0.7, 3, 1]))
        [0.22082609516300733, 0.7049003266226588]
    """
    # Checking input values
    max_species_number = 2
    max_inputs_per_species = 3
    assert len(crop_list) == max_species_number, "Only two species allowed"
    assert (
        len(crop_list[0]) == max_inputs_per_species
    ), "Only 3 inputs per species: ext_coeff, LAI, height"
    assert (
        len(crop_list[1]) == max_inputs_per_species
    ), "Only 3 inputs per species: ext_coeff, LAI, height"

    # Read inputs
    (
        extinction_coeff1,
        extinction_coeff2,
        height1,
        height2,
        leaf_area_index1,
        leaf_area_index2,
    ) = _read_rad_intercpt_wallace_inputs(crop_list=crop_list)

    transm_rad1 = math.exp(-extinction_coeff1 * leaf_area_index1)
    transm_rad2 = math.exp(-extinction_coeff2 * leaf_area_index2)

    # Dominant species rad interception
    rad_intercpt_dom1 = 1 - transm_rad1
    rad_intercpt_dom2 = 1 - transm_rad2
    height_fraction1 = height1 / (height1 + height2)
    height_fraction2 = height2 / (height1 + height2)

    # Suppressed species rad interception
    rad_intercpt_suppr1 = rad_intercpt_dom1 * transm_rad2
    rad_intercpt_suppr2 = rad_intercpt_dom2 * transm_rad1

    # Species rad interception
    rad_intercpt1 = rad_intercpt_suppr1 + height_fraction1 * (
        rad_intercpt_dom1 - rad_intercpt_suppr1
    )
    rad_intercpt2 = rad_intercpt_suppr2 + height_fraction2 * (
        rad_intercpt_dom2 - rad_intercpt_suppr2
    )
    return [rad_intercpt1, rad_intercpt2]


def _read_rad_intercpt_wallace_inputs(
    crop_list: tuple[list[float | int], list[float | int]],
):
    extinction_coeff1 = crop_list[0][0]
    extinction_coeff2 = crop_list[1][0]
    leaf_area_index1 = crop_list[0][1]
    leaf_area_index2 = crop_list[1][1]
    height1 = crop_list[0][2]
    height2 = crop_list[1][2]
    return (
        extinction_coeff1,
        extinction_coeff2,
        height1,
        height2,
        leaf_area_index1,
        leaf_area_index2,
    )


def rad_intercpt_apsim(
    crop_list: tuple[list[float | int], list[float | int]],
) -> NDArray:
    """Returns rad intercepted on each species.

    Args:
        crop_list: list of species characteristics:
        extinction_coeff: rad extinction coefficient
        leaf_area_index: leaf area index (m2/m2)

    References:
        Carberry, P.S., Adiku, S.G.K., McCown, R.L., Keating, B.A., 1996.Application of
         the APSIM cropping systems model to intercropping systems, in: Ito, C.,
         Johansen, C., Adu-Gyamfi, K., Katayama, K., Kumar-Rao, J.V.D.K., Rego,
         T.J. (Eds.), Dynamics of roots and nitrogen in cropping systems of the
         semi-arid tropics. Japan Int. Res. Centre Agric.Sci, pp. 637-648.

    Examples:
        >>> rad_intercpt_apsim(([0.5, 1],[0.7, 3]))
        array([0.17802431, 0.74770211])
    """
    # Variables init
    number_species = len(crop_list)
    transm_rad_temp = np.zeros(number_species)
    rad_intercpt_temp = np.zeros(number_species)
    ext_coeff_leaf_area_index_prod = np.zeros(number_species)
    cum_transm_rad = 1.0  # Cumulative fractional transmission
    rad_intercpt_list = np.zeros(number_species)

    # Temporarly rad interception
    for i in range(number_species):
        extinction_coeff = crop_list[i][0]
        leaf_area_index = crop_list[i][1]
        ext_coeff_leaf_area_index_prod[i] = extinction_coeff * leaf_area_index
        transm_rad_temp[i] = math.exp(-extinction_coeff * leaf_area_index)
        cum_transm_rad *= transm_rad_temp[i]
        rad_intercpt_temp[i] = 1 - transm_rad_temp[i]
    tot_rad_intercpt = 1 - cum_transm_rad

    # Actual rad interception of each crop weighted by k and LAI
    for i in range(number_species):
        if rad_intercpt_temp[i] > 0:
            rad_intercpt = (
                tot_rad_intercpt
                * ext_coeff_leaf_area_index_prod[i]
                / ext_coeff_leaf_area_index_prod.sum()
            )
        else:
            rad_intercpt = 0
        rad_intercpt_list[i] = rad_intercpt
    return rad_intercpt_list
