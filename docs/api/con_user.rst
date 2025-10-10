con_user
========

User constraint integration module.

Module Overview
---------------

The ``con_user`` module provides the interface for integrating user-defined
constraints into the trajectory optimization problem.

.. note::
   This module requires ``user_constraints.py`` in the working directory with
   ``equality_user()`` and ``inequality_user()`` functions defined.

Key Functions:
    * **equality_con_user**: Wrapper for user-defined equality constraints
    * **inequality_con_user**: Wrapper for user-defined inequality constraints
    * **equality_jac_con_user**: Jacobian of user equality constraints
    * **inequality_jac_con_user**: Jacobian of user inequality constraints

See ``example/user_constraints.py`` for a template.
