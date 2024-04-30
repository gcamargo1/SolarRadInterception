import numpy as np

from environ_biophysics import (
    height_weight_fact,
    opt_air_mass,
    rad_ext_coeff_black_beam,
    rad_ext_coeff_black_diff,
    solar_beam_fraction,
    solar_diffuse_fraction,
    rad_intercpt_sub_daily,
)


def test_height_weight_fact():
    assert height_weight_fact(0.75, 1.52, 0.56, 3) == 0.89
    assert height_weight_fact(1.5, 1.52, 0.56, 3) == 1.13


def test_opt_air_mass():
    assert opt_air_mass(100, 50) == 1.5357589603755304
    assert opt_air_mass(91.6, 30) == 1.0441319774485631


def test_rad_ext_coeff_black_beam():
    assert rad_ext_coeff_black_beam(0.087, 0) == 0.055576591963547875
    assert rad_ext_coeff_black_beam(0.087, 2) == 0.7254823957447912


def test_rad_ext_coeff_black_diff():
    assert rad_ext_coeff_black_diff(2, 0.1) == 0.9710451784887358
    assert rad_ext_coeff_black_diff(0, 0.1) == 0.9099461266386164


def test_solar_beam_fraction():
    assert solar_beam_fraction(101.3, 0, 0.75) == 0.75
    assert solar_beam_fraction(101.3, 50, 0.45) == 0.2892276326469122


def test_solar_diffuse_fraction():
    assert solar_diffuse_fraction(101.3, 0, 0.75) == 0.075
    assert solar_diffuse_fraction(101.3, 50, 0.45) == 0.10606799311188815


def test_rad_intercpt_sub_daily():
    result = rad_intercpt_sub_daily(
        0.75,
        101.3,
        0.8,
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
        ],
        0.5,
        2,
        np.linspace(0, 90, 19),
    )
    assert round(result[0][0], 4) == 0.0031
    assert round(result[0][1], 5) == 0.16775
