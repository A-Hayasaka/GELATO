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
Dynamics Module
===============

Module for defining rocket equations of motion.

This module implements rocket equations of motion in 6 degrees of freedom
(position, velocity, attitude), considering thrust, aerodynamic drag, and gravity.

Main Features:
    * Calculate time derivative of velocity (acceleration)
    * Equations of motion with aerodynamics
    * Equations of motion without aerodynamics (vacuum)
    * Time derivative of attitude using quaternions
    * Time derivative of mass (propellant consumption)

Functions:
    dynamics_velocity: Velocity equations of motion (with aerodynamics)
    dynamics_velocity_NoAir: Velocity equations of motion (without aerodynamics)
    dynamics_quaternion: Attitude equations of motion
    dynamics_mass: Mass equations of motion
"""

import numpy as np
from numpy.linalg import norm
from .coordinate import (
    quatrot,
    conj,
    quatmult,
    vel_eci2ecef,
    quat_nedg2eci,
    ecef2geodetic,
    ecef2eci,
    gravity,
)
from .USStandardAtmosphere import (
    airdensity_at,
    airpressure_at,
    geopotential_altitude,
    speed_of_sound,
)
from .utils import wind_ned


def dynamics_velocity(
    mass_e, pos_eci_e, vel_eci_e, quat_eci2body, t, param, wind_table, CA_table, units
):
    """Calculate time derivative of velocity (acceleration) with aerodynamics.
    
    Computes acceleration in ECI frame considering thrust, gravity, and
    aerodynamic forces (drag). Includes atmospheric density, pressure, wind,
    and Mach number effects.
    
    Args:
        mass_e (ndarray): Mass array (scaled) [kg]
        pos_eci_e (ndarray): Position in ECI frame (scaled) [m]
        vel_eci_e (ndarray): Velocity in ECI frame (scaled) [m/s]
        quat_eci2body (ndarray): Quaternion from ECI to body frame
        t (ndarray): Time array [s]
        param (ndarray): Vehicle parameters [thrust, Isp, area, etc.]
        wind_table (ndarray): Wind velocity lookup table
        CA_table (ndarray): Axial force coefficient vs Mach number
        units (ndarray): Scaling factors for mass, position, velocity
    
    Returns:
        ndarray: Acceleration in ECI frame (scaled) [m/s²]
    """

    mass = mass_e * units[0]
    pos_eci = pos_eci_e * units[1]
    vel_eci = vel_eci_e * units[2]
    acc_eci = np.zeros(vel_eci_e.shape)

    thrust_vac = param[0]
    air_area = param[2]
    nozzle_area = param[4]

    for i in range(len(mass)):
        pos_llh = ecef2geodetic(pos_eci[i, 0], pos_eci[i, 1], pos_eci[i, 2])
        altitude = geopotential_altitude(pos_llh[2])
        rho = airdensity_at(altitude)
        p = airpressure_at(altitude)

        vel_ecef = vel_eci2ecef(vel_eci[i], pos_eci[i], t[i])
        vel_wind_ned = wind_ned(altitude, wind_table)

        vel_wind_eci = quatrot(quat_nedg2eci(pos_eci[i], t[i]), vel_wind_ned)
        vel_air_eci = ecef2eci(vel_ecef, t[i]) - vel_wind_eci
        mach_number = norm(vel_air_eci) / speed_of_sound(altitude)

        ca = np.interp(mach_number, CA_table[:, 0], CA_table[:, 1])

        aeroforce_eci = 0.5 * rho * norm(vel_air_eci) * -vel_air_eci * air_area * ca

        thrust = thrust_vac - nozzle_area * p
        thrustdir_eci = quatrot(conj(quat_eci2body[i]), np.array([1.0, 0.0, 0.0]))
        thrust_eci = thrustdir_eci * thrust
        gravity_eci = gravity(pos_eci[i])

        acc_eci[i] = gravity_eci + (thrust_eci + aeroforce_eci) / mass[i]

    return acc_eci / units[2]


def dynamics_velocity_NoAir(mass_e, pos_eci_e, quat_eci2body, param, units):
    """Calculate time derivative of velocity (acceleration) without aerodynamics.
    
    Computes acceleration in ECI frame considering only thrust and gravity,
    without aerodynamic forces. Used for vacuum flight or when aerodynamics
    are negligible.
    
    Args:
        mass_e (ndarray): Mass array (scaled) [kg]
        pos_eci_e (ndarray): Position in ECI frame (scaled) [m]
        quat_eci2body (ndarray): Quaternion from ECI to body frame
        param (ndarray): Vehicle parameters [thrust, Isp, area, etc.]
        units (ndarray): Scaling factors for mass, position, velocity
    
    Returns:
        ndarray: Acceleration in ECI frame (scaled) [m/s²]
    """

    mass = mass_e * units[0]
    pos_eci = pos_eci_e * units[1]
    acc_eci = np.zeros(pos_eci_e.shape)

    thrust_vac = param[0]

    for i in range(len(mass)):

        thrust = thrust_vac
        thrustdir_eci = quatrot(conj(quat_eci2body[i]), np.array([1.0, 0.0, 0.0]))
        thrust_eci = thrustdir_eci * thrust
        gravity_eci = gravity(pos_eci[i])

        acc_eci[i] = gravity_eci + (thrust_eci) / mass[i]

    return acc_eci / units[2]


def dynamics_quaternion(quat_eci2body, u_e, unit_u):
    """Calculate time derivative of attitude quaternion.
    
    Computes the rate of change of quaternion based on angular velocity
    using the kinematic differential equation: dq/dt = 0.5 * q ⊗ ω
    
    Args:
        quat_eci2body (ndarray): Quaternion from ECI to body frame
        u_e (ndarray): Angular velocity control input (scaled) [deg/s]
        unit_u (float): Scaling factor for angular velocity
    
    Returns:
        ndarray: Time derivative of quaternion
    """

    u = u_e * unit_u

    d_quat = np.zeros(quat_eci2body.shape)
    for i in range(len(u)):
        omega_rps_body = np.deg2rad(np.array([0.0, u[i, 0], u[i, 1], u[i, 2]]))
        d_quat[i] = 0.5 * quatmult(quat_eci2body[i], omega_rps_body)

    return d_quat
