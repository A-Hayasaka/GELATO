Trajectory_Optimization
=======================

Main optimization program for solving launch trajectory problems.

.. note::
   This module requires ``user_constraints.py`` to be present in the working directory.
   See the ``example/`` directory for a template.

Module Overview
---------------

The ``Trajectory_Optimization`` module is the main entry point for solving trajectory
optimization problems using the Legendre-Gauss-Radau pseudospectral method.

Key Features:
    * Multi-phase trajectory optimization with LGR collocation
    * Flexible constraint formulation (equality and inequality)
    * Support for user-defined constraints
    * Automatic differentiation using finite differences
    * Interior point optimizer (IPOPT) interface

Usage:
    1. Prepare settings JSON file with vehicle parameters and constraints
    2. Create ``user_constraints.py`` with custom constraint functions
    3. Run optimization::

        python Trajectory_Optimization.py settings.json

    Results are saved to the output directory in CSV format.
