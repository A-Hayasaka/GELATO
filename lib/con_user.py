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
User Constraints Module
========================

Module for integrating user-defined constraint conditions.

This module incorporates user-defined equality and inequality constraints
described in user_constraints.py into optimization problems.

Main Features:
    * Integration of user-defined equality constraints
    * Integration of user-defined inequality constraints
    * Jacobian calculation using finite difference method
    * Interface provision for custom constraints

Functions:
    equality_user_wrapper: Wrapper for user-defined equality constraints
    inequality_user_wrapper: Wrapper for user-defined inequality constraints
"""

# constraints_u.py
# constraints about user conditions

from user_constraints import equality_user, inequality_user
from .jac_fd import jac_fd


def equality_jac_user(xdict, pdict, unitdict, condition):
    """Compute Jacobian of user-defined equality constraints using finite differences.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with problem parameters
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict or None: Jacobian matrices computed via finite differences, or None if no user constraints defined
    """
    if equality_user(xdict, pdict, unitdict, condition) is not None:
        return jac_fd(equality_user, xdict, pdict, unitdict, condition)


def inequality_jac_user(xdict, pdict, unitdict, condition):
    """Compute Jacobian of user-defined inequality constraints using finite differences.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with problem parameters
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict or None: Jacobian matrices computed via finite differences, or None if no user constraints defined
    """
    if inequality_user(xdict, pdict, unitdict, condition) is not None:
        return jac_fd(inequality_user, xdict, pdict, unitdict, condition)
