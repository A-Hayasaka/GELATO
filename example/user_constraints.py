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

import numpy as np
from lib.coordinate_c import orbital_elements
from lib.usercon_tools import get_value


def equality_user(xdict, pdict, unitdict, condition):
    """
    set additional  equality constraints.

    the return values of this function will be constrained to zero.

    it is strongly recommended to normalize the return values.

    the type of return value is float64 or numpy.ndarray(float64)

    """

    # get state value vector
    pos_IIP0 = get_value(xdict, pdict, unitdict, "IIP_END", "position")
    vel_IIP0 = get_value(xdict, pdict, unitdict, "IIP_END", "velocity")

    elem = orbital_elements(pos_IIP0, vel_IIP0)
    ha = (elem[0] * (1.0 - elem[1]) / 6378137.0) - 1.0  # height of apogee at the event IIP_END

    return ha


def equality_sparsity_user(pdict, condition):
    """
    set sparsity pattern of equality_user function.
    """

    return {
        "position": {"IIP_END": [0, 1, 2]},
        "velocity": {"IIP_END": [0, 1, 2]},
    }


def inequality_user(xdict, pdict, unitdict, condition):
    """
    set additional inequality constraints.

    the return values of this function will be constrained to positive or zero.

    it is strongly recommended to normalize the return values.

    the type of return value is float64 or numpy.ndarray(float64)

    """

    return None


def inequality_sparsity_user(pdict, condition):
    """
    set sparsity pattern of inequality_user function.
    """

    return None
