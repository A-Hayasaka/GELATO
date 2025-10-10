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
Section Parameters Module
==========================

Module for managing parameters of each phase (section).

This module centrally manages pseudospectral method parameters (collocation points,
differentiation matrices, etc.) for each section in multi-phase optimization problems.

Main Features:
    * Manage LGR collocation points for each section
    * Manage differentiation matrices for each section
    * Convert indices between sections
    * Structure variable arrays

Classes:
    PSparams: Management class for pseudospectral method parameters
        - Holds number of sections and number of nodes in each section
        - Pre-calculates LGR collocation points and differentiation matrices
        - Provides index conversion methods
"""

import numpy as np
from .PSfunctions import differentiation_matrix_LGR, nodes_LGR


class PSparams:
    """Pseudospectral method parameters management class.
    
    Manages collocation points, differentiation matrices, and index mappings
    for multi-section pseudospectral optimization problems using Legendre-
    Gauss-Radau (LGR) collocation.
    
    Attributes:
        _num_sections (int): Number of trajectory sections
        _num_nodes (list): Number of LGR nodes per section
        _tau (list): LGR collocation points for each section
        _D (list): Differentiation matrices for each section
        _index_start_u (list): Starting indices for control variables
        _N (int): Total number of collocation nodes
    """
    
    def __init__(self, num_nodes: "list[int]"):
        """Initialize PSparams with node counts for each section.
        
        Args:
            num_nodes (list[int]): Number of LGR nodes for each section
        """
        self._num_sections = len(num_nodes)
        self._num_nodes = [n for n in num_nodes]
        self._tau = [nodes_LGR(n) for n in num_nodes]
        self._D = [differentiation_matrix_LGR(n) for n in num_nodes]
        self._index_start_u = [
            sum(self._num_nodes[:i]) for i in range(self._num_sections)
        ]
        self._N = sum(self._num_nodes)

    def tau(self, i):
        """Get LGR collocation points for section i.
        
        Args:
            i (int): Section index
        
        Returns:
            ndarray: LGR collocation points τ ∈ [-1, 1]
        """
        if i < 0 or i >= self._num_sections:
            raise ValueError("Index out of range")
        return self._tau[i]

    def D(self, i):
        """Get differentiation matrix for section i.
        
        Args:
            i (int): Section index
        
        Returns:
            ndarray: LGR differentiation matrix D
        """
        if i < 0 or i >= self._num_sections:
            raise ValueError("Index out of range")
        return self._D[i]

    def index_start_u(self, i):
        """Get starting index for control variables in section i.
        
        Args:
            i (int): Section index
        
        Returns:
            int: Starting index in flattened control array
        """
        return self._index_start_u[i]

    def index_end_u(self, i):
        """Get ending index for control variables in section i.
        
        Args:
            i (int): Section index
        
        Returns:
            int: Ending index in flattened control array
        """
        return self.index_start_u(i) + self._num_nodes[i]

    def index_start_x(self, i):
        """Get starting index for state variables in section i.
        
        Args:
            i (int): Section index
        
        Returns:
            int: Starting index in flattened state array
        """
        return self._index_start_u[i] + i

    def index_end_x(self, i):
        """Get ending index for state variables in section i.
        
        Args:
            i (int): Section index
        
        Returns:
            int: Ending index in flattened state array
        """
        return self.index_start_x(i) + self._num_nodes[i] + 1

    def num_u(self):
        """Get total number of control points across all sections.
        
        Returns:
            int: Total number of LGR collocation nodes
        """
        return self._N

    def num_x(self):
        """Get total number of state points across all sections.
        
        Returns:
            int: Total number of state nodes (including section boundaries)
        """
        return self._N + self._num_sections

    def num_sections(self):
        """Get number of sections.
        
        Returns:
            int: Number of trajectory sections
        """
        return self._num_sections

    def nodes(self, i):
        """Get number of LGR nodes in section i.
        
        Args:
            i (int): Section index
        
        Returns:
            int: Number of collocation nodes
        """
        if i < 0 or i >= self._num_sections:
            raise ValueError("Index out of range")
        return self._num_nodes[i]

    def time_nodes(self, i, to, tf):
        """Get physical time values at collocation points for section i.
        
        Converts normalized LGR points τ ∈ [-1, 1] to physical time
        t ∈ [to, tf] using affine transformation.
        
        Args:
            i (int): Section index
            to (float): Initial time of section [s]
            tf (float): Final time of section [s]
        
        Returns:
            ndarray: Time values at collocation points [s]
        """
        t = np.zeros(self._num_nodes[i] + 1)
        t[0] = to
        t[1:] = self.tau(i) * (tf - to) / 2 + (tf + to) / 2
        return t

    def get_index(self, section):
        """get set of index for a given section.

        Args:
            section (int): section number

        Returns:
            ua (int): start index for u
            ub (int): end index for u
            xa (int): start index for x
            xb (int): end index for x
            n (int): number of nodes in the section
        """

        ua = self._index_start_u[section]
        n = self._num_nodes[section]
        ub = ua + n
        xa = ua + section
        xb = xa + n + 1

        return ua, ub, xa, xb, n

    def __getitem__(self, i):
        """Get parameters dictionary for section i (backward compatibility).
        
        Args:
            i (int): Section index
        
        Returns:
            dict: Dictionary with keys 'index_start', 'nodes', 'D', 'tau'
        """
        if i < 0 or i >= self._num_sections:
            raise ValueError("Index out of range")
        return {
            "index_start": self._index_start_u[i],
            "nodes": self._num_nodes[i],
            "D": self._D[i],
            "tau": self._tau[i],
        }
