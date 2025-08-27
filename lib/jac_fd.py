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
def jac_fd(con, xdict, pdict, unitdict, condition, sparsity="ALL"):
    """
    Calculate jacobian by finite-difference method(forward difference).

    Note that this function is slow, because this function do not use sparse matrix.

    Args:
        con(function) : object function
        xdict : variable arg for con
        pdict : parameter arg for con
        unitdict : unit arg for con
        conditiondict : condition arg for con
        sparsity(dict) : sparsity pattern for jacobian (optional)

    Returns:
        jac(dict(ndarray)) : dict of jacobian matrix

    """

    jac = {}
    dx = pdict["dx"]
    g_base = con(xdict, pdict, unitdict, condition)
    if hasattr(g_base, "__len__"):
        nRows = len(g_base)
    else:
        nRows = 1

    @profile
    def subroutine(key, i):
        xdict[key][i] += dx
        g_p = con(xdict, pdict, unitdict, condition)
        xdict[key][i] -= dx
        print(i, g_p - g_base)
        if nRows == 1:
            jac[key]["coo"][0].append(0)
            jac[key]["coo"][1].append(i)
            jac[key]["coo"][2].append((g_p - g_base) / dx)
        else:
            jac[key]["coo"][0].extend(range(nRows))
            jac[key]["coo"][1].extend([i] * nRows)
            jac[key]["coo"][2].extend((g_p - g_base) / dx)

    for key, val in xdict.items():
        jac[key] = {
            "coo": [[], [], []],
            "shape": (nRows, val.size)
        }
        if sparsity == "ALL":
            for i in range(val.size):
                subroutine(key, i)

        elif key in sparsity:
            if key == "t":
                for i in range(val.size):
                    subroutine(key, i)
            else:
                for event_index in range(pdict["ps_params"].num_sections()):
                    xa = pdict["ps_params"].index_start_x(event_index)
                    xb = pdict["ps_params"].index_end_x(event_index)
                    print(event_index, xa, xb)
                    for i in range(xa * 3, xb * 3):
                        subroutine(key, i)

        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac
