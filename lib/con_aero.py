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
Aerodynamic Constraints Module
===============================

Module for defining aerodynamic constraint conditions.

This module calculates aerodynamic constraints (dynamic pressure, heating rate,
angle of attack, etc.) during rocket flight and provides them as constraint
conditions for optimization problems.

Main Features:
    * Calculate and constrain dynamic pressure
    * Calculate q-α (dynamic pressure × angle of attack) constraints
    * Calculate heating rate
    * Angle of attack constraints
    * Equality constraints calculation
    * Inequality constraints calculation

Constraint Functions:
    dynamic_pressure_dimless: Non-dimensionalized dynamic pressure
    angle_of_attack_dimless: Non-dimensionalized angle of attack
    q_alpha_dimless: Non-dimensionalized q-α product
    heating_rate_dimless: Non-dimensionalized heating rate
"""

# constraints_c.py
# constraints about aerodynamic conditions

import numpy as np
from .utils_c import (
    dynamic_pressure_pa,
    angle_of_attack_all_rad,
    q_alpha_pa_rad,
    dynamic_pressure_array_pa,
    angle_of_attack_all_array_rad,
    q_alpha_array_pa_rad,
)
from .coordinate_c import *


def dynamic_pressure_dimless(pos_eci_e, vel_eci_e, t_e, wind, units):
    """Calculate normalized dynamic pressure.
    
    Computes dynamic pressure q = 0.5 * ρ * v² normalized by maximum
    allowable value for constraint evaluation.
    
    Args:
        pos_eci_e (ndarray): Position in ECI frame (scaled)
        vel_eci_e (ndarray): Velocity in ECI frame (scaled)
        t_e (float): Time (scaled)
        wind (ndarray): Wind velocity table
        units (ndarray): [pos_scale, vel_scale, time_scale, q_max]
    
    Returns:
        float: Normalized dynamic pressure (dimensionless)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return dynamic_pressure_pa(pos_eci, vel_eci, t, wind) / units[3]


def dynamic_pressure_array_dimless(pos_eci_e, vel_eci_e, t_e, wind, units):
    """
    Calculate dimensionless dynamic pressure array for trajectory points.
    
    Computes dynamic pressure (0.5 * ρ * v²) at each state point and normalizes
    by the characteristic pressure unit for constraint evaluation.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized position array in ECI frame (n, 3)
        vel_eci_e (numpy.ndarray): Normalized velocity array in ECI frame (n, 3)
        t_e (numpy.ndarray): Normalized time array (n,)
        wind (numpy.ndarray): Wind profile data table
        units (tuple): Unit scaling factors (pos_unit, vel_unit, time_unit, pressure_unit)
    
    Returns:
        numpy.ndarray: Dimensionless dynamic pressure at each point (n,)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return dynamic_pressure_array_pa(pos_eci, vel_eci, t, wind) / units[3]


def angle_of_attack_all_dimless(pos_eci_e, vel_eci_e, quat, t_e, wind, units):
    """Calculate normalized total angle of attack.
    
    Computes total angle of attack (magnitude combining pitch and yaw)
    normalized by maximum allowable value for constraint evaluation.
    
    Args:
        pos_eci_e (ndarray): Position in ECI frame (scaled)
        vel_eci_e (ndarray): Velocity in ECI frame (scaled)
        quat (ndarray): Attitude quaternion
        t_e (float): Time (scaled)
        wind (ndarray): Wind velocity table
        units (ndarray): [pos_scale, vel_scale, time_scale, alpha_max]
    
    Returns:
        float: Normalized total angle of attack (dimensionless)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return angle_of_attack_all_rad(pos_eci, vel_eci, quat, t, wind) / units[3]


def angle_of_attack_all_array_dimless(pos_eci_e, vel_eci_e, quat, t_e, wind, units):
    """
    Calculate dimensionless total angle of attack array for trajectory points.
    
    Computes total angle of attack (magnitude combining pitch and yaw components)
    at each state point and normalizes by maximum allowable angle for constraint evaluation.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized position array in ECI frame (n, 3)
        vel_eci_e (numpy.ndarray): Normalized velocity array in ECI frame (n, 3)
        quat (numpy.ndarray): Attitude quaternion array (n, 4)
        t_e (numpy.ndarray): Normalized time array (n,)
        wind (numpy.ndarray): Wind profile data table
        units (tuple): Unit scaling factors (pos_unit, vel_unit, time_unit, alpha_max)
    
    Returns:
        numpy.ndarray: Dimensionless angle of attack at each point (n,)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return angle_of_attack_all_array_rad(pos_eci, vel_eci, quat, t, wind) / units[3]


def q_alpha_dimless(pos_eci_e, vel_eci_e, quat, t_e, wind, units):
    """Calculate normalized Q-alpha product.
    
    Computes Q-alpha (dynamic pressure × angle of attack) normalized
    by maximum allowable value. Critical constraint for structural loads.
    
    Args:
        pos_eci_e (ndarray): Position in ECI frame (scaled)
        vel_eci_e (ndarray): Velocity in ECI frame (scaled)
        quat (ndarray): Attitude quaternion
        t_e (float): Time (scaled)
        wind (ndarray): Wind velocity table
        units (ndarray): [pos_scale, vel_scale, time_scale, q_alpha_max]
    
    Returns:
        float: Normalized Q-alpha (dimensionless)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return q_alpha_pa_rad(pos_eci, vel_eci, quat, t, wind) / units[3]


def q_alpha_array_dimless(pos_eci_e, vel_eci_e, quat, t_e, wind, units):
    """
    Calculate dimensionless Q-alpha array for trajectory points.
    
    Computes Q-alpha (dynamic pressure × angle of attack) at each state point
    and normalizes by maximum allowable value. This parameter is critical for
    assessing aerodynamic structural loads during flight.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized position array in ECI frame (n, 3)
        vel_eci_e (numpy.ndarray): Normalized velocity array in ECI frame (n, 3)
        quat (numpy.ndarray): Attitude quaternion array (n, 4)
        t_e (numpy.ndarray): Normalized time array (n,)
        wind (numpy.ndarray): Wind profile data table
        units (tuple): Unit scaling factors (pos_unit, vel_unit, time_unit, q_alpha_max)
    
    Returns:
        numpy.ndarray: Dimensionless Q-alpha at each point (n,)
    """
    pos_eci = pos_eci_e * units[0]
    vel_eci = vel_eci_e * units[1]
    t = t_e * units[2]
    return q_alpha_array_pa_rad(pos_eci, vel_eci, quat, t, wind) / units[3]


def inequality_max_alpha(xdict, pdict, unitdict, condition):
    """
    Inequality constraint for maximum angle of attack.
    
    Enforces that angle of attack remains below specified maximum value
    throughout flight or during specified trajectory sections. Constraint
    is formulated as: 1 - α/α_max ≥ 0.
    
    Args:
        xdict (dict): State variables (position, velocity, quaternion, time)
        pdict (dict): Problem parameters (sections, pseudospectral params, wind)
        unitdict (dict): Unit scaling factors
        condition (dict): Constraint conditions including 'AOA_max' specifications
    
    Returns:
        numpy.ndarray: Constraint values (should be ≥ 0)
    """

    con = []

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    t = xdict["t"]

    num_sections = pdict["num_sections"]

    wind = pdict["wind_table"]

    for i in range(num_sections - 1):
        section_name = pdict["params"][i]["name"]

        # angle of attack
        if section_name in condition["AOA_max"]:

            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]
            quat_i_ = quat_[xa:xb]

            to = t[i]
            tf = t[i + 1]

            aoa_max = condition["AOA_max"][section_name]["value"] * np.pi / 180.0
            units[3] = aoa_max
            if condition["AOA_max"][section_name]["range"] == "all":
                t_i_ = pdict["ps_params"].time_nodes(i, to, tf)
                con.append(
                    1.0
                    - angle_of_attack_all_array_dimless(
                        pos_i_, vel_i_, quat_i_, t_i_, wind, units
                    )
                )
            elif condition["AOA_max"][section_name]["range"] == "initial":
                con.append(
                    1.0
                    - angle_of_attack_all_dimless(
                        pos_i_[0], vel_i_[0], quat_i_[0], to, wind, units
                    )
                )

    if len(con) == 0:
        return None
    else:
        return np.concatenate(con, axis=None)


def inequality_max_q(xdict, pdict, unitdict, condition):
    """Inequality constraint enforcing maximum dynamic pressure limits.
    
    Constrains dynamic pressure Q = 0.5 * ρ * V² to remain below specified maximum values
    during specified flight phases.
    
    Args:
        xdict (dict): Dictionary containing state variables ('position', 'velocity', 't')
        pdict (dict): Dictionary with problem parameters including wind table
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary with 'dynamic_pressure_max' specifications
    
    Returns:
        ndarray or None: Constraint values (>= 0 when satisfied), or None if no constraints
    """

    con = []

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)

    t = xdict["t"]

    num_sections = pdict["num_sections"]

    wind = pdict["wind_table"]

    for i in range(num_sections - 1):
        section_name = pdict["params"][i]["name"]
        # max-Q
        if section_name in condition["dynamic_pressure_max"]:

            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]
            to = t[i]
            tf = t[i + 1]
            q_max = condition["dynamic_pressure_max"][section_name]["value"]
            units[3] = q_max
            if condition["dynamic_pressure_max"][section_name]["range"] == "all":
                t_i_ = pdict["ps_params"].time_nodes(i, to, tf)
                con.append(
                    1.0
                    - dynamic_pressure_array_dimless(pos_i_, vel_i_, t_i_, wind, units)
                )
            elif condition["dynamic_pressure_max"][section_name]["range"] == "initial":
                con.append(
                    1.0
                    - dynamic_pressure_dimless(pos_i_[0], vel_i_[0], to, wind, units)
                )

    if len(con) == 0:
        return None
    else:
        return np.concatenate(con, axis=None)


def inequality_max_qalpha(xdict, pdict, unitdict, condition):
    """Inequality constraint enforcing maximum Q-alpha limit.
    
    Q-alpha is the product of dynamic pressure and angle of attack, representing
    combined aerothermal and structural loading. This constraint ensures the vehicle
    stays within safe aerodynamic heating and load limits during flight.
    
    The constraint can be applied to all collocation points or only the initial point
    of specified trajectory sections.
    
    Args:
        xdict (dict): Dictionary with 'position', 'velocity', 'quaternion', 't' state variables
        pdict (dict): Dictionary with 'num_sections', 'params', 'ps_params', 'wind_table'
        unitdict (dict): Dictionary with 'position', 'velocity', 't' unit scalings
        condition (dict): Configuration with 'Q_alpha_max' limits per section
            - 'value': Maximum Q-alpha limit in [Pa·rad]
            - 'range': 'all' (all points) or 'initial' (section start only)
    
    Returns:
        numpy.ndarray or None: Constraint values (>= 0 when satisfied), or None if no constraints
    """

    con = []

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    t = xdict["t"]

    num_sections = pdict["num_sections"]

    wind = pdict["wind_table"]

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]
        # max-Qalpha
        if section_name in condition["Q_alpha_max"]:

            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]
            quat_i_ = quat_[xa:xb]
            to = t[i]
            tf = t[i + 1]

            qalpha_max = condition["Q_alpha_max"][section_name]["value"] * np.pi / 180.0
            units[3] = qalpha_max
            if condition["Q_alpha_max"][section_name]["range"] == "all":
                t_i_ = pdict["ps_params"].time_nodes(i, to, tf)
                con.append(
                    1.0
                    - q_alpha_array_dimless(pos_i_, vel_i_, quat_i_, t_i_, wind, units)
                )
            elif condition["Q_alpha_max"][section_name]["range"] == "initial":
                con.append(
                    1.0
                    - q_alpha_dimless(pos_i_[0], vel_i_[0], quat_i_[0], to, wind, units)
                )

    if len(con) == 0:
        return None
    else:
        return np.concatenate(con, axis=None)


def inequality_length_max_alpha(xdict, pdict, unitdict, condition):
    """Calculate the number of maximum angle-of-attack constraints.
    
    Determines constraint vector length for Jacobian matrix sizing.
    
    Returns:
        int: Total number of AOA constraint equations
    """
    res = 0

    num_sections = pdict["num_sections"]

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]
        # max-Qalpha
        if section_name in condition["AOA_max"]:

            if condition["AOA_max"][section_name]["range"] == "all":
                res += pdict["ps_params"].nodes(i) + 1
            elif condition["AOA_max"][section_name]["range"] == "initial":
                res += 1

    return res


def inequality_length_max_q(xdict, pdict, unitdict, condition):
    """Calculate the number of maximum dynamic pressure constraints.
    
    Determines constraint vector length for Jacobian matrix sizing.
    
    Returns:
        int: Total number of dynamic pressure constraint equations
    """
    res = 0

    num_sections = pdict["num_sections"]

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]
        # max-Qalpha
        if section_name in condition["dynamic_pressure_max"]:

            if condition["dynamic_pressure_max"][section_name]["range"] == "all":
                res += pdict["ps_params"].nodes(i) + 1
            elif condition["dynamic_pressure_max"][section_name]["range"] == "initial":
                res += 1

    return res


def inequality_length_max_qalpha(xdict, pdict, unitdict, condition):
    """Calculate the number of maximum Q-alpha constraints.
    
    Determines constraint vector length for Jacobian matrix sizing.
    Q-alpha = dynamic pressure × angle of attack.
    
    Returns:
        int: Total number of Q-alpha constraint equations
    """
    res = 0

    num_sections = pdict["num_sections"]

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]
        # max-Qalpha
        if section_name in condition["Q_alpha_max"]:

            if condition["Q_alpha_max"][section_name]["range"] == "all":
                res += pdict["ps_params"].nodes(i) + 1
            elif condition["Q_alpha_max"][section_name]["range"] == "initial":
                res += 1

    return res


def angle_of_attack_all_gradient_dimless(
    pos_eci_e, vel_eci_e, quat, to, tf, wind, units, time_nodes, dx, n
):
    """Compute gradient of angle of attack with respect to state variables using finite differences.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized ECI position array
        vel_eci_e (numpy.ndarray): Normalized ECI velocity array
        quat (numpy.ndarray): Quaternion array (ECI to body frame)
        to (float): Normalized initial time for section
        tf (float): Normalized final time for section
        wind (object): Wind table for atmospheric calculations
        units (list): Unit scaling factors [position, velocity, time, dimensionless]
        time_nodes (callable): Function to generate collocation time nodes
        dx (float): Finite difference step size
        n (int): Number of collocation points
    
    Returns:
        dict: Gradient dictionary with keys 'position', 'velocity', 'quaternion', 'to', 'tf'
    """
    ki = range(n)
    grad = {
        "position": np.zeros((n, 3)),
        "velocity": np.zeros((n, 3)),
        "quaternion": np.zeros((n, 4)),
        "to": np.zeros(n),
        "tf": np.zeros(n),
    }

    pos_ki = pos_eci_e[ki]
    vel_ki = vel_eci_e[ki]
    quat_ki = quat[ki]

    t_nodes = time_nodes(to, tf)
    t_ki = t_nodes[ki]

    f_c = angle_of_attack_all_array_dimless(pos_ki, vel_ki, quat_ki, t_ki, wind, units)

    for j in range(3):
        pos_ki[:, j] += dx
        f_p = angle_of_attack_all_array_dimless(
            pos_ki, vel_ki, quat_ki, t_ki, wind, units
        )
        pos_ki[:, j] -= dx
        grad["position"][:, j] = (f_p - f_c) / dx

    for j in range(3):
        vel_ki[:, j] += dx
        f_p = angle_of_attack_all_array_dimless(
            pos_ki, vel_ki, quat_ki, t_ki, wind, units
        )
        vel_ki[:, j] -= dx
        grad["velocity"][:, j] = (f_p - f_c) / dx

    for j in range(4):
        quat_ki[:, j] += dx
        f_p = angle_of_attack_all_array_dimless(
            pos_ki, vel_ki, quat_ki, t_ki, wind, units
        )
        quat_ki[:, j] -= dx
        grad["quaternion"][:, j] = (f_p - f_c) / dx

    to_p = to + dx
    t_nodes_p1 = time_nodes(to_p, tf)
    f_p = angle_of_attack_all_array_dimless(
        pos_ki, vel_ki, quat_ki, t_nodes_p1[ki], wind, units
    )
    grad["to"] = (f_p - f_c) / dx

    tf_p = tf + dx
    t_nodes_p2 = time_nodes(to, tf_p)
    f_p = angle_of_attack_all_array_dimless(
        pos_ki, vel_ki, quat_ki, t_nodes_p2[ki], wind, units
    )
    grad["tf"] = (f_p - f_c) / dx

    return grad


def inequality_jac_max_alpha(xdict, pdict, unitdict, condition):
    """Jacobian of the maximum angle of attack inequality constraint.
    
    Computes derivatives of inequality_max_alpha with respect to state variables.
    
    Args:
        xdict (dict): Dictionary with 'position', 'velocity', 'quaternion', 't'
        pdict (dict): Dictionary with 'dx' (FD step), 'wind_table', 'num_sections', 'params'
        unitdict (dict): Dictionary with 'position', 'velocity', 't' unit scalings
        condition (dict): Configuration dictionary with alpha constraint settings
    
    Returns:
        dict or None: Jacobian in COO sparse format with keys 'data', 'row', 'col', or None if no constraints
    """

    jac = {}
    dx = pdict["dx"]

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    t = xdict["t"]

    wind = pdict["wind_table"]
    num_sections = pdict["num_sections"]

    nRow = inequality_length_max_alpha(xdict, pdict, unitdict, condition)
    if nRow == 0:
        return None

    jac["position"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 4)}
    jac["t"] = {"coo": [[], [], []], "shape": (nRow, num_sections + 1)}

    iRow = 0

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]

        # angle of attack
        if section_name in condition["AOA_max"]:
            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]
            quat_i_ = quat_[xa:xb]

            to = t[i]
            tf = t[i + 1]

            aoa_max = condition["AOA_max"][section_name]["value"] * np.pi / 180.0
            units[3] = aoa_max

            if condition["AOA_max"][section_name]["range"] == "all":
                nk = n + 1
            elif condition["AOA_max"][section_name]["range"] == "initial":
                nk = 1

            def time_nodes(t1, t2):
                """Generate collocation time nodes for current section."""
                return pdict["ps_params"].time_nodes(i, t1, t2)

            dfdx = angle_of_attack_all_gradient_dimless(
                pos_i_, vel_i_, quat_i_, to, tf, wind, units, time_nodes, dx, nk
            )

            for j in range(3):
                jac["position"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["position"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["position"]["coo"][2].extend(-dfdx["position"].ravel(order="F"))

            for j in range(3):
                jac["velocity"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["velocity"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["velocity"]["coo"][2].extend(-dfdx["velocity"].ravel(order="F"))

            for j in range(4):
                jac["quaternion"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["quaternion"]["coo"][1].extend(
                    range(xa * 4 + j, (xa + nk) * 4 + j, 4)
                )
            jac["quaternion"]["coo"][2].extend(-dfdx["quaternion"].ravel(order="F"))

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i] * nk)
            jac["t"]["coo"][2].extend(-dfdx["to"])

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i + 1] * nk)
            jac["t"]["coo"][2].extend(-dfdx["tf"])

            iRow += nk

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac


def dynamic_pressure_gradient_dimless(
    pos_eci_e, vel_eci_e, to, tf, wind, units, time_nodes, dx, n
):
    """Compute gradient of dynamic pressure with respect to state variables using finite differences.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized ECI position array
        vel_eci_e (numpy.ndarray): Normalized ECI velocity array
        to (float): Normalized initial time for section
        tf (float): Normalized final time for section
        wind (object): Wind table for atmospheric calculations
        units (list): Unit scaling factors [position, velocity, time, dimensionless]
        time_nodes (callable): Function to generate collocation time nodes
        dx (float): Finite difference step size
        n (int): Number of collocation points
    
    Returns:
        dict: Gradient dictionary with keys 'position', 'velocity', 'to', 'tf'
    """
    ki = range(n)
    grad = {
        "position": np.zeros((n, 3)),
        "velocity": np.zeros((n, 3)),
        "to": np.zeros(n),
        "tf": np.zeros(n),
    }

    pos_ki = pos_eci_e[ki]
    vel_ki = vel_eci_e[ki]

    t_nodes = time_nodes(to, tf)
    t_ki = t_nodes[ki]

    f_c = dynamic_pressure_array_dimless(pos_ki, vel_ki, t_ki, wind, units)

    for j in range(3):
        pos_ki[:, j] += dx
        f_p = dynamic_pressure_array_dimless(pos_ki, vel_ki, t_ki, wind, units)
        pos_ki[:, j] -= dx
        grad["position"][:, j] = (f_p - f_c) / dx

    for j in range(3):
        vel_ki[:, j] += dx
        f_p = dynamic_pressure_array_dimless(pos_ki, vel_ki, t_ki, wind, units)
        vel_ki[:, j] -= dx
        grad["velocity"][:, j] = (f_p - f_c) / dx

    to_p = to + dx
    t_nodes_p1 = time_nodes(to_p, tf)
    f_p = dynamic_pressure_array_dimless(pos_ki, vel_ki, t_nodes_p1[ki], wind, units)
    grad["to"] = (f_p - f_c) / dx

    tf_p = tf + dx
    t_nodes_p2 = time_nodes(to, tf_p)
    f_p = dynamic_pressure_array_dimless(pos_ki, vel_ki, t_nodes_p2[ki], wind, units)
    grad["tf"] = (f_p - f_c) / dx

    return grad


def inequality_jac_max_q(xdict, pdict, unitdict, condition):
    """Jacobian of the maximum dynamic pressure inequality constraint.
    
    Computes derivatives of inequality_max_q with respect to state variables.
    
    Args:
        xdict (dict): Dictionary with 'position', 'velocity', 't' state variables
        pdict (dict): Dictionary with 'dx' (FD step), 'wind_table', 'num_sections', 'params'
        unitdict (dict): Dictionary with 'position', 'velocity', 't' unit scalings
        condition (dict): Configuration dictionary with dynamic pressure constraint settings
    
    Returns:
        dict or None: Jacobian in COO sparse format with keys 'data', 'row', 'col', or None if no constraints
    """

    jac = {}
    dx = pdict["dx"]

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)

    t = xdict["t"]

    wind = pdict["wind_table"]
    num_sections = pdict["num_sections"]

    nRow = inequality_length_max_q(xdict, pdict, unitdict, condition)
    if nRow == 0:
        return None

    jac["position"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 4)}
    jac["t"] = {"coo": [[], [], []], "shape": (nRow, num_sections + 1)}

    iRow = 0

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]

        # angle of attack
        if section_name in condition["dynamic_pressure_max"]:
            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]

            to = t[i]
            tf = t[i + 1]

            q_max = condition["dynamic_pressure_max"][section_name]["value"]
            units[3] = q_max

            if condition["dynamic_pressure_max"][section_name]["range"] == "all":
                nk = n + 1
            elif condition["dynamic_pressure_max"][section_name]["range"] == "initial":
                nk = 1

            def time_nodes(t1, t2):
                """Generate collocation time nodes for current section."""
                return pdict["ps_params"].time_nodes(i, t1, t2)

            dfdx = dynamic_pressure_gradient_dimless(
                pos_i_, vel_i_, to, tf, wind, units, time_nodes, dx, nk
            )

            for j in range(3):
                jac["position"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["position"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["position"]["coo"][2].extend(-dfdx["position"].ravel(order="F"))

            for j in range(3):
                jac["velocity"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["velocity"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["velocity"]["coo"][2].extend(-dfdx["velocity"].ravel(order="F"))

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i] * nk)
            jac["t"]["coo"][2].extend(-dfdx["to"])

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i + 1] * nk)
            jac["t"]["coo"][2].extend(-dfdx["tf"])

            iRow += nk

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac


def q_alpha_gradient_dimless(
    pos_eci_e, vel_eci_e, quat, to, tf, wind, units, time_nodes, dx, n
):
    """Compute gradient of Q-alpha (dynamic pressure × angle of attack) using finite differences.
    
    Args:
        pos_eci_e (numpy.ndarray): Normalized ECI position array
        vel_eci_e (numpy.ndarray): Normalized ECI velocity array
        quat (numpy.ndarray): Quaternion array (ECI to body frame)
        to (float): Normalized initial time for section
        tf (float): Normalized final time for section
        wind (object): Wind table for atmospheric calculations
        units (list): Unit scaling factors [position, velocity, time, dimensionless]
        time_nodes (callable): Function to generate collocation time nodes
        dx (float): Finite difference step size
        n (int): Number of collocation points
    
    Returns:
        dict: Gradient dictionary with keys 'position', 'velocity', 'quaternion', 'to', 'tf'
    """
    ki = range(n)
    grad = {
        "position": np.zeros((n, 3)),
        "velocity": np.zeros((n, 3)),
        "quaternion": np.zeros((n, 4)),
        "to": np.zeros(n),
        "tf": np.zeros(n),
    }

    pos_ki = pos_eci_e[ki]
    vel_ki = vel_eci_e[ki]
    quat_ki = quat[ki]

    t_nodes = time_nodes(to, tf)
    t_ki = t_nodes[ki]

    f_c = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_ki, wind, units)

    for j in range(3):
        pos_ki[:, j] += dx
        f_p = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_ki, wind, units)
        pos_ki[:, j] -= dx
        grad["position"][:, j] = (f_p - f_c) / dx

    for j in range(3):
        vel_ki[:, j] += dx
        f_p = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_ki, wind, units)
        vel_ki[:, j] -= dx
        grad["velocity"][:, j] = (f_p - f_c) / dx

    for j in range(4):
        quat_ki[:, j] += dx
        f_p = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_ki, wind, units)
        quat_ki[:, j] -= dx
        grad["quaternion"][:, j] = (f_p - f_c) / dx

    to_p = to + dx
    t_nodes_p1 = time_nodes(to_p, tf)
    f_p = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_nodes_p1[ki], wind, units)
    grad["to"] = (f_p - f_c) / dx

    tf_p = tf + dx
    t_nodes_p2 = time_nodes(to, tf_p)
    f_p = q_alpha_array_dimless(pos_ki, vel_ki, quat_ki, t_nodes_p2[ki], wind, units)
    grad["tf"] = (f_p - f_c) / dx

    return grad


def inequality_jac_max_qalpha(xdict, pdict, unitdict, condition):
    """Jacobian of the maximum Q-alpha (dynamic pressure × angle of attack) inequality constraint.
    
    Computes derivatives of inequality_max_qalpha with respect to state variables.
    
    Args:
        xdict (dict): Dictionary with 'position', 'velocity', 'quaternion', 't' state variables
        pdict (dict): Dictionary with 'dx' (FD step), 'wind_table', 'num_sections', 'params'
        unitdict (dict): Dictionary with 'position', 'velocity', 't' unit scalings
        condition (dict): Configuration dictionary with Q-alpha constraint settings
    
    Returns:
        dict or None: Jacobian in COO sparse format with keys 'data', 'row', 'col', or None if no constraints
    """

    jac = {}
    dx = pdict["dx"]

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]
    unit_t = unitdict["t"]
    units = [unit_pos, unit_vel, unit_t, 1.0]

    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    t = xdict["t"]

    wind = pdict["wind_table"]
    num_sections = pdict["num_sections"]

    nRow = inequality_length_max_qalpha(xdict, pdict, unitdict, condition)
    if nRow == 0:
        return None

    jac["position"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 4)}
    jac["t"] = {"coo": [[], [], []], "shape": (nRow, num_sections + 1)}

    iRow = 0

    for i in range(num_sections - 1):

        section_name = pdict["params"][i]["name"]

        # angle of attack
        if section_name in condition["Q_alpha_max"]:
            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            pos_i_ = pos_[xa:xb]
            vel_i_ = vel_[xa:xb]
            quat_i_ = quat_[xa:xb]

            to = t[i]
            tf = t[i + 1]

            qalpha_max = condition["Q_alpha_max"][section_name]["value"] * np.pi / 180.0
            units[3] = qalpha_max

            if condition["Q_alpha_max"][section_name]["range"] == "all":
                nk = n + 1
            elif condition["Q_alpha_max"][section_name]["range"] == "initial":
                nk = 1

            def time_nodes(t1, t2):
                """Generate collocation time nodes for current section."""
                return pdict["ps_params"].time_nodes(i, t1, t2)

            dfdx = q_alpha_gradient_dimless(
                pos_i_, vel_i_, quat_i_, to, tf, wind, units, time_nodes, dx, nk
            )
            for j in range(3):
                jac["position"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["position"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["position"]["coo"][2].extend(-dfdx["position"].ravel(order="F"))

            for j in range(3):
                jac["velocity"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["velocity"]["coo"][1].extend(
                    range(xa * 3 + j, (xa + nk) * 3 + j, 3)
                )
            jac["velocity"]["coo"][2].extend(-dfdx["velocity"].ravel(order="F"))

            for j in range(4):
                jac["quaternion"]["coo"][0].extend(range(iRow, iRow + nk))
                jac["quaternion"]["coo"][1].extend(
                    range(xa * 4 + j, (xa + nk) * 4 + j, 4)
                )
            jac["quaternion"]["coo"][2].extend(-dfdx["quaternion"].ravel(order="F"))

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i] * nk)
            jac["t"]["coo"][2].extend(-dfdx["to"])

            jac["t"]["coo"][0].extend(range(iRow, iRow + nk))
            jac["t"]["coo"][1].extend([i + 1] * nk)
            jac["t"]["coo"][2].extend(-dfdx["tf"])

            iRow += nk

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac
