import numpy as np

from environ_biophysics import (
    get_height_weight_factor,
    get_optical_air_mass,
    get_solar_radiation_extinction_coeff_black_beam,
    get_solar_radiation_extinction_coeff_black_diff,
    get_solar_beam_fraction,
    get_solar_diffuse_fraction,
    get_solar_radiation_interception_sub_daily,
)


def test_get_height_weight_factor():
    assert get_height_weight_factor(0.75, 1.52, 0.56, 3) == 0.89
    assert get_height_weight_factor(1.5, 1.52, 0.56, 3) == 1.13


def test_get_optical_air_mass():
    assert get_optical_air_mass(100, 50) == 1.5357589603755304
    assert get_optical_air_mass(91.6, 30) == 1.0441319774485631


def test_get_solar_radiation_extinction_coeff_black_beam():
    assert (
        get_solar_radiation_extinction_coeff_black_beam(0.087, 0)
        == 0.055576591963547875
    )
    assert (
        get_solar_radiation_extinction_coeff_black_beam(0.087, 2) == 0.7254823957447912
    )


def test_get_solar_radiation_extinction_coeff_black_diff():
    assert get_solar_radiation_extinction_coeff_black_diff(2, 0.1) == 0.9710451784887358
    assert get_solar_radiation_extinction_coeff_black_diff(0, 0.1) == 0.9099461266386164


def test_get_solar_beam_fraction():
    assert get_solar_beam_fraction(101.3, 0, 0.75) == 0.75
    assert get_solar_beam_fraction(101.3, 50, 0.45) == 0.2892276326469122


def test_get_solar_diffuse_fraction():
    assert get_solar_diffuse_fraction(101.3, 0, 0.75) == 0.075
    assert get_solar_diffuse_fraction(101.3, 50, 0.45) == 0.10606799311188815


def test_get_solar_radiation_interception_sub_daily():
    result = get_solar_radiation_interception_sub_daily(
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
