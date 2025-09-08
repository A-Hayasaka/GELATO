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

@profile
def jac_fd(con, xdict, pdict, unitdict, condition, keycol_nonzero="ALL"):
    """
    Calculate jacobian by finite-difference method(forward difference).

    Note that this function is slow, because this function do not use sparse matrix.

    Args:
        con(function) : object function
        xdict : variable arg for con
        pdict : parameter arg for con
        unitdict : unit arg for con
        conditiondict : condition arg for con
        keycol_nonzero(dict or "ALL") : sparsity pattern of jacobian matrix

    Returns:
        jac(dict(ndarray)) : dict of jacobian matrix

    """

    jac = {}
    dx = pdict["dx"]
    num_sections = pdict["num_sections"]
    g_base = con(xdict, pdict, unitdict, condition)
    if hasattr(g_base, "__len__"):
        nRows = len(g_base)
    else:
        nRows = 1
    jac["mass"] = {"coo": [[], [], []], "shape": (nRows, pdict["M"])}
    jac["position"] = {"coo": [[], [], []], "shape": (nRows, pdict["M"] * 3)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRows, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRows, pdict["M"] * 4)}
    jac["u"] = {"coo": [[], [], []], "shape": (nRows, pdict["N"] * 3)}
    jac["t"] = {"coo": [[], [], []], "shape": (nRows, num_sections + 1)}

    for key, val in xdict.items():
        for i in range(val.size):
            xdict[key][i] += dx
            g_p = con(xdict, pdict, unitdict, condition)
            jac[key]["coo"][0].extend(list(range(nRows)))
            jac[key]["coo"][1].extend([i] * nRows)
            jac[key]["coo"][2].extend((g_p - g_base) / dx)
            xdict[key][i] -= dx

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype=np.int32)
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype=np.int32)
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype=np.float64)

    return jac
