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
Utilities Module
================

Module providing general-purpose helper functions.

This module provides auxiliary functions for calculating various physical
quantities and coordinate transformations.

Main Features:
    * Calculate dynamic pressure
    * Calculate angle of attack
    * Calculate heating rate
    * Interpolate wind velocity
    * Helper functions for vector operations

Functions:
    dynamic_pressure_pa: Calculate dynamic pressure [Pa]
    angle_of_attack_all_rad: Calculate total angle of attack [rad]
    q_alpha_pa_rad: Calculate q-α product
    heating_rate: Calculate heating rate [W/m^2]
    wind_ned: Get wind velocity in NED coordinate system
"""

from math import sin, cos, asin, atan2, sqrt, radians, degrees
import numpy as np
from numpy.linalg import norm
from .USStandardAtmosphere import *
from .coordinate import *


def haversine(lon1, lat1, lon2, lat2, r):
    """Calculate great-circle distance using haversine formula (DEPRECATED).
    
    Computes spherical distance between two points. For accurate geodetic
    distance on WGS84 ellipsoid, use Vincenty formula instead.
    
    Args:
        lon1 (float64): Longitude of start point [deg]
        lat1 (float64): Latitude of start point [deg]
        lon2 (float64): Longitude of end point [deg]
        lat2 (float64): Latitude of end point [deg]
        r (float64): Sphere radius [m]
    
    Returns:
        float64: Great-circle distance [m]
    
    Note:
        DEPRECATED - Use Vincenty formula for WGS84 ellipsoid accuracy
    """

    # convert decimal degrees to radians
    lon1 = radians(lon1)
    lat1 = radians(lat1)
    lon2 = radians(lon2)
    lat2 = radians(lat2)

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    return 2 * r * asin(sqrt(a))


def distance_on_earth(x_km, y_km, z_km, lon0, lat0, time):
    """Calculate surface distance using haversine formula (DEPRECATED).
    
    Computes distance from reference point to ECEF position using spherical
    approximation. Use Vincenty formula for accurate WGS84 ellipsoid distance.
    
    Args:
        x_km (float64): X-coordinate of target in ECEF frame [km]
        y_km (float64): Y-coordinate of target in ECEF frame [km]
        z_km (float64): Z-coordinate of target in ECEF frame [km]
        lon0 (float64): Longitude of reference point [deg]
        lat0 (float64): Latitude of reference point [deg]
        time (float64): Time [s]
    
    Returns:
        float64: Surface distance [km]
    
    Note:
        DEPRECATED - Use Vincenty formula for WGS84 ellipsoid accuracy
    """
    radius_km = 6378.142
    angular_velocity_rps = 0.729211586e-4

    r_km = norm(np.array([x_km, y_km, z_km]))
    latitude = asin(z_km / r_km)
    longitude = atan2(y_km, x_km) - angular_velocity_rps * time

    return haversine(lon0, lat0, degrees(longitude), degrees(latitude), radius_km)


def wind_ned(altitude_m, wind_data):
    """Get wind velocity in NED frame by linear interpolation.
    
    Args:
        altitude_m (float): Altitude [m]
        wind_data (ndarray): Wind table with columns [altitude, north, east]
    
    Returns:
        ndarray: Wind velocity [north, east, down] in NED frame [m/s]
    """
    wind = np.zeros(3)
    wind[0] = np.interp(altitude_m, wind_data[:, 0], wind_data[:, 1])
    wind[1] = np.interp(altitude_m, wind_data[:, 0], wind_data[:, 2])
    wind[2] = 0.0
    return wind


def angle_of_attack_all_rad(pos_eci, vel_eci, quat, t, wind):
    """Calculate total angle of attack.
    
    Computes angle between velocity vector (relative to air) and vehicle
    thrust axis. Accounts for wind effects.
    
    Args:
        pos_eci (ndarray): Position in ECI frame [m]
        vel_eci (ndarray): Inertial velocity in ECI frame [m/s]
        quat (ndarray): Quaternion from ECI to body frame
        t (float64): Time [s]
        wind (ndarray): Wind velocity table
    
    Returns:
        float64: Total angle of attack [rad]
    """

    thrust_dir_eci = quatrot(conj(quat), np.array([1.0, 0.0, 0.0]))

    pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2])
    altitude_m = geopotential_altitude(pos_llh[2])

    vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t)
    vel_wind_ned = wind_ned(altitude_m, wind)

    vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned)
    vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci

    c_alpha = normalize(vel_air_eci).dot(normalize(thrust_dir_eci))

    if c_alpha >= 1.0 or norm(vel_air_eci) < 0.001:
        return 0.0
    else:
        return acos(c_alpha)


def angle_of_attack_all_array_rad(pos_eci, vel_eci, quat, t, wind):
    """Calculate total angle of attack for array of states.
    
    Vectorized version of angle_of_attack_all_rad.
    
    Args:
        pos_eci (ndarray): Position array in ECI frame [m]
        vel_eci (ndarray): Velocity array in ECI frame [m/s]
        quat (ndarray): Quaternion array from ECI to body frame
        t (ndarray): Time array [s]
        wind (ndarray): Wind velocity table
    
    Returns:
        ndarray: Total angle of attack array [rad]
    """
    alpha = np.zeros(pos_eci.shape[0])
    for i in range(pos_eci.shape[0]):
        alpha[i] = angle_of_attack_all_rad(pos_eci[i], vel_eci[i], quat[i], t[i], wind)
    return alpha


def angle_of_attack_ab_rad(pos_eci, vel_eci, quat, t, wind):
    """Calculate pitch and yaw angles of attack separately.
    
    Decomposes total angle of attack into pitch (alpha) and yaw (beta)
    components in body frame.
    
    Args:
        pos_eci (ndarray): Position in ECI frame [m]
        vel_eci (ndarray): Inertial velocity in ECI frame [m/s]
        quat (ndarray): Quaternion from ECI to body frame
        t (float64): Time [s]
        wind (ndarray): Wind velocity table
    
    Returns:
        ndarray: [alpha, beta] - pitch and yaw angles of attack [rad]
    """

    pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2])
    altitude_m = geopotential_altitude(pos_llh[2])

    vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t)
    vel_wind_ned = wind_ned(altitude_m, wind)

    vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned)
    vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci

    vel_air_body = quatrot(quat, vel_air_eci)

    if vel_air_body[0] < 0.001:
        return np.zeros(2)
    else:
        alpha_z = atan2(vel_air_body[2], vel_air_body[0])
        alpha_y = atan2(vel_air_body[1], vel_air_body[0])
        return np.array((alpha_z, alpha_y))


def dynamic_pressure_pa(pos_eci, vel_eci, t, wind):
    """Calculates dynamic pressure.
    Args:
        pos_eci (ndarray) : position in ECI frame [m]
        vel_eci (ndarray) : inertial velocity in ECI frame [m/s]
        t (float64) : time [s]
        wind (ndarray) : wind table
    Returns:
        float64 : dynamic pressure [Pa]
    """

    pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2])
    altitude_m = geopotential_altitude(pos_llh[2])
    rho = airdensity_at(altitude_m)

    vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t)
    vel_wind_ned = wind_ned(altitude_m, wind)
    vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned)
    vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci

    return 0.5 * vel_air_eci.dot(vel_air_eci) * rho


def dynamic_pressure_array_pa(pos_eci, vel_eci, t, wind):
    """Calculate dynamic pressure for multiple trajectory points.
    
    Vectorized version of dynamic_pressure_pa for processing trajectory arrays.
    
    Args:
        pos_eci (ndarray): Array of ECI positions (N, 3) [m]
        vel_eci (ndarray): Array of ECI velocities (N, 3) [m/s]
        t (ndarray): Array of times (N,) [s]
        wind (callable): Wind model function
    
    Returns:
        ndarray: Dynamic pressure values (N,) [Pa]
    """
    q = np.zeros(pos_eci.shape[0])
    for i in range(pos_eci.shape[0]):
        q[i] = dynamic_pressure_pa(pos_eci[i], vel_eci[i], t[i], wind)
    return q


def q_alpha_pa_rad(pos_eci, vel_eci, quat, t, wind):
    """Calculate Q-alpha (dynamic pressure × angle of attack).
    
    Args:
        pos_eci (ndarray): ECI position (3,) [m]
        vel_eci (ndarray): ECI velocity (3,) [m/s]
        quat (ndarray): Attitude quaternion (4,)
        t (float): Time [s]
        wind (callable): Wind model function
    
    Returns:
        float: Q-alpha [Pa·rad]
    """
    alpha = angle_of_attack_all_rad(pos_eci, vel_eci, quat, t, wind)
    q = dynamic_pressure_pa(pos_eci, vel_eci, t, wind)
    return q * alpha


def q_alpha_array_pa_rad(pos_eci, vel_eci, quat, t, wind):
    """Calculate Q-alpha for multiple trajectory points.
    
    Vectorized version of q_alpha_pa_rad for processing trajectory arrays.
    
    Args:
        pos_eci (ndarray): Array of ECI positions (N, 3) [m]
        vel_eci (ndarray): Array of ECI velocities (N, 3) [m/s]
        quat (ndarray): Array of quaternions (N, 4)
        t (ndarray): Array of times (N,) [s]
        wind (callable): Wind model function
    
    Returns:
        ndarray: Q-alpha values (N,) [Pa·rad]
    """
    qa = np.zeros(pos_eci.shape[0])
    for i in range(pos_eci.shape[0]):
        qa[i] = q_alpha_pa_rad(pos_eci[i], vel_eci[i], quat[i], t[i], wind)
    return qa
