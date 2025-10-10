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
Pseudospectral Method Functions Module
=======================================

Module providing implementation of Legendre-Gauss-Radau pseudospectral method.

This module provides the mathematical foundation for pseudospectral methods
to discretize trajectory optimization problems. Uses LGR (Legendre-Gauss-Radau)
collocation points to achieve high-precision numerical calculations.

Main Features:
    * Calculate LGR collocation points (nodes)
    * Calculate Legendre polynomials
    * Generate differentiation matrices
    * Lagrange interpolation polynomials
    * Calculate integration weights

References:
    * Benson, D. A. (2005). A Gauss Pseudospectral Transcription for Optimal Control
    * Garg et al. (2009). An overview of three pseudospectral methods

Functions:
    LegendreFunction: Calculate Legendre polynomials
    nodes_LGR: Calculate LGR collocation points
    differentiation_matrix_LGR: Generate LGR differentiation matrix
    lagrange: Lagrange interpolation polynomial
"""

from scipy import special
import numpy as np


def LegendreFunction(x, n):
    """Evaluate Legendre polynomial of the first kind Pn(x).

    Args:
        x (float64): Argument of the Legendre function, typically ∈ [-1, 1]
        n (int): Degree of the Legendre polynomial

    Returns:
        float64: Value of Pn(x)
    """

    Legendre, Derivative = special.lpn(n, x)
    return Legendre[-1]


def lagrange(tn, k, t):
    """Evaluate kth Lagrange basis polynomial.
    
    Computes Lk(t) = ∏(i≠k) (t - ti)/(tk - ti) for interpolation.

    Args:
        tn (ndarray): Collocation points (nodes)
        k (int): Index of the basis polynomial
        t (float64): Evaluation point

    Returns:
        float64: Value of Lk(t)
    """
    L = 1.0
    N = len(tn)
    for i in range(N):
        if i != k:
            L = L * (t - tn[i]) / (tn[k] - tn[i])
    return L


def lagrangeD(tn, k, t):
    """Evaluate derivative of kth Lagrange basis polynomial.
    
    Computes dLk/dt at point t using analytic differentiation.

    Args:
        tn (ndarray): Collocation points (nodes)
        k (int): Index of the basis polynomial
        t (float64): Evaluation point

    Returns:
        float64: Value of dLk/dt at t
    """
    N = len(tn)
    den = 1.0
    for i in range(N):
        if i != k:
            den = den * (tn[k] - tn[i])
    num = 0.0
    for j in range(N):
        num_j = 1.0
        if j != k:
            for i in range(N):
                if i != k and i != j:
                    num_j = num_j * (t - tn[i])
            num = num + num_j
    return num / den


def nodes_LGL(n):
    """Compute Legendre-Gauss-Lobatto (LGL) collocation points.
    
    Returns n points including boundary points at -1 and +1.
    
    Args:
        n (int): Number of LGL points
    
    Returns:
        ndarray: LGL nodes in [-1, 1]
    """
    roots, weight = special.j_roots(n - 2, 1, 1)
    nodes = np.hstack((-1, roots, 1))
    return nodes


def weight_LGL(n):
    """Compute Legendre-Gauss-Lobatto (LGL) quadrature weights.
    
    Args:
        n (int): Number of LGL points
    
    Returns:
        ndarray: LGL quadrature weights
    """
    nodes = nodes_LGL(n)
    w = np.zeros(0)
    for i in range(n):
        w = np.append(w, 2 / (n * (n - 1) * LegendreFunction(nodes[i], n - 1) ** 2))
    return w


def differentiation_matrix_LGL(n):
    """Compute Legendre-Gauss-Lobatto (LGL) differentiation matrix.
    
    Generates matrix D such that df/dτ ≈ D·f for interpolating polynomial.
    
    Args:
        n (int): Number of LGL points
    
    Returns:
        ndarray: n×n differentiation matrix
    """
    tau = nodes_LGL(n)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = (
                    LegendreFunction(tau[i], n - 1)
                    / LegendreFunction(tau[j], n - 1)
                    / (tau[i] - tau[j])
                )
            elif i == j and i == 0:
                D[i, j] = -n * (n - 1) * 0.25
            elif i == j and i == n - 1:
                D[i, j] = n * (n - 1) * 0.25
            else:
                D[i, j] = 0.0
    return D


def nodes_LG(n):
    """Compute Legendre-Gauss (LG) collocation points.
    
    Args:
        n (int): Number of LG points
    
    Returns:
        ndarray: LG nodes in (-1, 1) (excludes boundaries)
    """
    return special.roots_legendre(n)[0]


def weight_LG(n):
    """Compute Legendre-Gauss (LG) quadrature weights.
    
    Args:
        n (int): Number of LG points
    
    Returns:
        ndarray: LG quadrature weights
    """
    return special.roots_legendre(n)[1]


def differentiation_matrix_LG(n):
    """Compute Legendre-Gauss (LG) differentiation matrix.
    
    Generates matrix D such that df/dτ ≈ D·f for interpolating polynomial.
    
    Args:
        n (int): Number of LG points
    
    Returns:
        ndarray: n×(n+1) differentiation matrix
    """
    tk_lg, _ = special.roots_legendre(n)
    tk_lg = np.hstack((-1.0, tk_lg))
    D = np.zeros((n, n + 1))
    for k in range(1, n + 1):
        for i in range(n + 1):
            D[k - 1, i] = lagrangeD(tk_lg, i, tk_lg[k])
    return D


def nodes_LGR(n, reverse=True):
    """Legendre-Gauss-Radau(LGR) points.

    Args:
        n (int) : number of degrees. (n >= 2)
        reverse (boolean) : type of LGR points. The return value
        includes -1 when reverse is false and it includes +1 when
        reverse is true.

    Returns:
        ndarray: LGR points.

    """

    roots, weight = special.j_roots(n - 1, 0, 1)
    nodes = np.hstack((-1, roots))
    if reverse:
        return np.sort(-nodes)
    else:
        return nodes


def weight_LGR(n):
    """Compute Legendre-Gauss-Radau (LGR) quadrature weights.
    
    Args:
        n (int): Number of LGR points (n >= 2)
    
    Returns:
        ndarray: LGR quadrature weights
    """
    nodes = nodes_LGR(n)
    w = np.zeros(0)
    for i in range(n):
        w = np.append(
            w, (1 - nodes[i]) / (n * n * LegendreFunction(nodes[i], n - 1) ** 2)
        )
    return w


def differentiation_matrix_LGR(n, reverse=True):
    """Compute Legendre-Gauss-Radau (LGR) differentiation matrix.
    
    Generates matrix D such that df/dτ ≈ D·f at LGR collocation points,
    where f includes values at both interior nodes and one boundary.

    Args:
        n (int): Number of LGR points (n >= 2)
        reverse (boolean): LGR type. If True, includes boundary at +1;
                          if False, includes boundary at -1

    Returns:
        ndarray: n×(n+1) differentiation matrix

    """

    tk_lgr = nodes_LGR(n, reverse)
    if reverse:
        tk_lgr = np.hstack((-1.0, tk_lgr))
    else:
        tk_lgr = np.hstack((tk_lgr, 1.0))
    D = np.zeros((n, n + 1))
    for k in range(n):
        for i in range(n + 1):
            if reverse:
                D[k, i] = lagrangeD(tk_lgr, i, tk_lgr[k + 1])
            else:
                D[k, i] = lagrangeD(tk_lgr, i, tk_lgr[k])
    return D
