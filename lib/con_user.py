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

# constraints_u.py
# constraints about user conditions

import user_constraints as uc
from .jac_fd import jac_fd


def equality_user(xdict, pdict, unitdict, condition):
    """User-defined equality constraint."""
    if not hasattr(uc, "equality_user"):
        return None
    else:
        return uc.equality_user(xdict, pdict, unitdict, condition)


def equality_sparsity_user(pdict, condition):
    """Sparsity pattern of user-defined equality constraint."""
    if not hasattr(uc, "equality_sparsity_user"):
        return "ALL"
    else:
        return uc.equality_sparsity_user(pdict, condition)


def equality_jac_user(xdict, pdict, unitdict, condition):
    """Jacobian of user-defined equality constraint."""
    if not hasattr(uc, "equality_user"):
        return None
    elif uc.equality_user(xdict, pdict, unitdict, condition) is not None:
        keycol_nonzero = equality_sparsity_user(pdict, condition)
        return jac_fd(
            uc.equality_user, xdict, pdict, unitdict, condition, keycol_nonzero
        )


def inequality_user(xdict, pdict, unitdict, condition):
    """User-defined inequality constraint."""
    if not hasattr(uc, "inequality_user"):
        return None
    else:
        return uc.inequality_user(xdict, pdict, unitdict, condition)


def inequality_sparsity_user(pdict, condition):
    """Sparsity pattern of user-defined inequality constraint."""
    if not hasattr(uc, "inequality_sparsity_user"):
        return "ALL"
    else:
        return uc.inequality_sparsity_user(pdict, condition)


def inequality_jac_user(xdict, pdict, unitdict, condition):
    """Jacobian of user-defined inequality constraint."""
    if not hasattr(uc, "inequality_user"):
        return None
    if uc.inequality_user(xdict, pdict, unitdict, condition) is not None:
        keycol_nonzero = inequality_sparsity_user(pdict, condition)
        return jac_fd(
            uc.inequality_user, xdict, pdict, unitdict, condition, keycol_nonzero
        )
