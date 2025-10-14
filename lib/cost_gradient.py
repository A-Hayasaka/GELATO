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
Cost Gradient Module
====================

Module for defining objective functions and their gradients.

This module calculates the objective function (cost function) for trajectory
optimization problems and its gradient vector.

Main Optimization Modes:
    * Payload: Maximize payload mass (maximize initial mass)
    * Time: Minimize arrival time (maximize residual propellant)

Functions:
    cost_6DoF: Calculate objective function
    cost_jac: Calculate gradient vector of objective function
"""

import numpy as np


def cost_6DoF(xdict, condition):
    """
    Objective (cost) function for trajectory optimization.
    
    Defines the optimization objective based on the selected mode:
    - 'Payload' mode: Maximize initial mass (negative cost to minimize)
    - Other modes: Minimize final time (maximize remaining propellant)
    
    Args:
        xdict (dict): Dictionary containing state variables:
            - 'mass': Mass array (dimensionless)
            - 't': Time array (dimensionless)
        condition (dict): Configuration with 'OptimizationMode' specification
    
    Returns:
        float: Scalar cost value to be minimized by the optimizer.
    """
    if condition["OptimizationMode"] == "Payload":
        return -xdict["mass"][0]  # Maximize initial mass (dimensionless)
    else:
        return xdict["t"][-1]  # Minimize arrival time (= maximize remaining propellant)


def cost_jac(xdict, condition):
    """
    Compute gradient (Jacobian) of the objective function.
    
    Calculates the partial derivatives of the cost function with respect to 
    optimization variables (mass or time depending on mode).
    
    Args:
        xdict (dict): Dictionary containing state variables
        condition (dict): Configuration with 'OptimizationMode' specification
    
    Returns:
        dict: Gradient dictionary with keys:
            - 'mass': Gradient w.r.t. mass (Payload mode)
            - 't': Gradient w.r.t. time (other modes)
    """

    jac = {}
    if condition["OptimizationMode"] == "Payload":
        jac["mass"] = np.zeros(xdict["mass"].size)
        jac["mass"][0] = -1.0
    else:
        jac["t"] = np.zeros(xdict["t"].size)
        jac["t"][-1] = 1.0
    return jac
