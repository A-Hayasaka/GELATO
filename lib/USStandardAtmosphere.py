#
# The MIT License
#
# Copyright (c) 2022 Interstellar Technologies Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

"""
U.S. Standard Atmosphere Module
================================

Module implementing atmospheric model based on U.S. Standard Atmosphere 1976.

This module calculates atmospheric physical quantities (density, pressure,
temperature, speed of sound, etc.) as a function of altitude.
Used for rocket aerodynamic calculations.

Main Features:
    * Calculate atmospheric density vs. altitude
    * Calculate atmospheric pressure vs. altitude
    * Calculate temperature vs. altitude
    * Calculate speed of sound
    * Convert geopotential altitude

References:
    U.S. Standard Atmosphere, 1976
    NASA-TM-X-74335

Functions:
    geopotential_altitude: Calculate geopotential altitude
    airdensity_at: Calculate atmospheric density
    airpressure_at: Calculate atmospheric pressure
    temperature_at: Calculate temperature
    speed_of_sound: Calculate speed of sound
"""

import math
from math import sqrt
import numpy as np

# U.S. Standard Atmosphere 1976


def geopotential_altitude(z):
    """Calculate geopotential altitude from geometric altitude.
    
    Converts geometric altitude to geopotential altitude using the relationship
    defined in U.S. Standard Atmosphere 1976, Equation 18. For altitudes above
    86 km, returns the geometric altitude unchanged.
    
    Args:
        z (float64): Geometric altitude [m]
    
    Returns:
        float64: Geopotential altitude [m] (z <= 86000 m), or geometric altitude (z > 86000 m)
    
    Reference:
        U.S. Standard Atmosphere, 1976, NASA-TM-X-74335
    """

    r0 = 6356766
    g0 = 9.80665

    if z < 86000:
        return 1.0 * (r0 * z) / (r0 + z)
    else:
        return z


def us_standard_atmosphere_params_at(altitude_m):
    """Get atmospheric parameters at reference levels for given altitude.
    
    Returns the atmospheric parameters for the appropriate atmospheric layer
    based on the U.S. Standard Atmosphere 1976 model. Temperature gradients
    for layers from 91-110 km and above 120 km are approximate values.
    
    Args:
        altitude_m (float64): Geopotential altitude (Z <= 86000 m) or
                             geometric altitude (Z > 86000 m) [m]
    
    Returns:
        ndarray: Atmospheric parameters [6 elements]
            [0]: Reference altitude [m]
            [1]: Temperature gradient (lapse rate) [K/m]
            [2]: Reference temperature [K]
            [3]: Reference pressure [Pa]
            [4]: Specific gas constant [J/(kg·K)]
            [5]: Gravity acceleration [m/s²]
    """

    GRAVITY_ACC_CONST = 9.80665
    Rstar = 8314.32

    #  b <  7:[Hb, Lmb, Tmb, Pb, R]
    #  b >= 7:[Zb, Lb,  Tb,  Pb, R]
    PARAMS = [
        [0.0, -0.0065, 288.15, 101325.0, Rstar / 28.9644],
        [11000.0, 0.0, 216.65, 22632.0, Rstar / 28.9644],
        [20000.0, 0.001, 216.65, 5474.9, Rstar / 28.9644],
        [32000.0, 0.0028, 228.65, 868.02, Rstar / 28.9644],
        [47000.0, 0.0, 270.65, 110.91, Rstar / 28.9644],
        [51000.0, -0.0028, 270.65, 66.939, Rstar / 28.9644],
        [71000.0, -0.002, 214.65, 3.9564, Rstar / 28.9644],
        [86000.0, 0.0, 186.8673, 0.37338, Rstar / 28.9522],
        [91000.0, 0.0025, 186.8673, 0.15381, Rstar / 28.89],  # Lb is dummy
        [110000.0, 0.012, 240.0, 7.1042e-3, Rstar / 27.27],
        [120000.0, 0.012, 360.0, 2.5382e-3, Rstar / 26.20],  # Lb is dummy
    ]

    k = 0
    for i in range(len(PARAMS)):
        if altitude_m >= PARAMS[i][0]:
            k = i

    return np.append(np.array(PARAMS[k]), GRAVITY_ACC_CONST)


def airtemperature_at(altitude_m):
    """Calculate air temperature at given altitude.
    
    Computes atmospheric temperature using different models for different
    altitude regimes as defined in U.S. Standard Atmosphere 1976, Section 1.2.5.
    Uses linear gradient below 91 km, molecular temperature model 91-110 km,
    and exponential model above 120 km.
    
    Args:
        altitude_m (float64): Geopotential altitude (Z <= 86000 m) or
                             geometric altitude (Z > 86000 m) [m]
    
    Returns:
        float64: Air temperature [K]
    
    Reference:
        U.S. Standard Atmosphere, 1976, Section 1.2.5
    """

    air_params = us_standard_atmosphere_params_at(altitude_m)
    HAL = air_params[0]
    LR = air_params[1]
    T0 = air_params[2]
    P0 = air_params[3]
    R = air_params[4]
    GRAVITY = air_params[5]

    if altitude_m <= 91000:
        return T0 + LR * (altitude_m - HAL)
    elif altitude_m <= 110000:
        Tc = 263.1905
        A = -76.3232
        a = -19942.9
        return Tc + A * sqrt(1.0 - (altitude_m - 91000) ** 2 / a**2)
    elif altitude_m <= 120000:
        return T0 + LR * (altitude_m - HAL)
    else:
        Tinf = 1000.0
        r0 = 6356766
        xi = (altitude_m - HAL) * (r0 + HAL) / (r0 + altitude_m)
        return Tinf - (Tinf - T0) * np.exp(-0.01875 * 1e-3 * xi)


def airpressure_at(altitude_m):
    """Calculate air pressure at given altitude.
    
    Computes atmospheric pressure using the hydrostatic equation as defined
    in U.S. Standard Atmosphere 1976, Section 1.3.1. Uses exact formulas for
    isothermal and gradient layers below 86 km, and approximate fitting for
    altitudes above 86 km.
    
    Args:
        altitude_m (float64): Geopotential altitude (Z <= 86000 m) or
                             geometric altitude (Z > 86000 m) [m]
    
    Returns:
        float64: Air pressure [Pa]
    
    Reference:
        U.S. Standard Atmosphere, 1976, Section 1.3.1
    """

    air_params = us_standard_atmosphere_params_at(altitude_m)
    HAL = air_params[0]
    LR = air_params[1]
    T0 = air_params[2]
    P0 = air_params[3]
    R = air_params[4]
    GRAVITY = air_params[5]

    air_temperature = airtemperature_at(altitude_m)

    if math.fabs(LR) > 1.0e-10:
        air_pressure = P0 * ((T0 + LR * (altitude_m - HAL)) / T0) ** (GRAVITY / -LR / R)
    else:
        air_pressure = P0 * math.exp(GRAVITY / R * (HAL - altitude_m) / T0)

    return air_pressure


def airdensity_at(altitude_m):
    """Calculate air density at given altitude.
    
    Computes atmospheric mass density using the ideal gas law and hydrostatic
    equation as defined in U.S. Standard Atmosphere 1976, Section 1.3.4.
    Values above 86 km are approximations using fitted models.
    
    Args:
        altitude_m (float64): Geopotential altitude (Z <= 86000 m) or
                             geometric altitude (Z > 86000 m) [m]
    
    Returns:
        float64: Air mass density [kg/m³]
    
    Reference:
        U.S. Standard Atmosphere, 1976, Section 1.3.4
    """

    air_params = us_standard_atmosphere_params_at(altitude_m)
    HAL = air_params[0]
    LR = air_params[1]
    T0 = air_params[2]
    P0 = air_params[3]
    R = air_params[4]
    GRAVITY = air_params[5]

    air_temperature = airtemperature_at(altitude_m)

    if math.fabs(LR) > 1.0e-10:
        air_pressure = P0 * ((T0 + LR * (altitude_m - HAL)) / T0) ** (GRAVITY / -LR / R)
    else:
        air_pressure = P0 * math.exp(GRAVITY / R * (HAL - altitude_m) / T0)

    air_density = air_pressure / R / air_temperature
    return air_density


def speed_of_sound(altitude_m):
    """Calculate speed of sound at given altitude.
    
    Computes the speed of sound in air using temperature as defined in
    U.S. Standard Atmosphere 1976, Section 1.3.10. Uses the formula
    a = sqrt(γ * R * T) where γ is the ratio of specific heats.
    
    Args:
        altitude_m (float64): Geopotential altitude (Z <= 86000 m) or
                             geometric altitude (Z > 86000 m) [m]
    
    Returns:
        float64: Speed of sound [m/s]
    
    Reference:
        U.S. Standard Atmosphere, 1976, Section 1.3.10
    """

    air_params = us_standard_atmosphere_params_at(altitude_m)
    HAL = air_params[0]
    LR = air_params[1]
    T0 = air_params[2]
    P0 = air_params[3]
    R = air_params[4]
    GRAVITY = air_params[5]

    air_temperature = airtemperature_at(altitude_m)

    return sqrt(1.4 * R * air_temperature)
