Installation
============

Requirements
------------

GELATO requires the following software:

* Python 3.8 or higher
* NumPy
* SciPy
* Matplotlib
* C++ compiler (for building PyBind11 modules)

Installation Steps
------------------

1. Clone the Repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/istellartech/GELATO.git
   cd GELATO

2. Install Dependencies
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install numpy scipy matplotlib

3. Build C++ Extension Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GELATO uses C++ modules bound with PyBind11 for performance optimization.
Build them with the following command:

.. code-block:: bash

   make

This will build the following C++ extension modules:

* ``coordinate_c``
* ``dynamics_c``
* ``IIP_c``
* ``USStandardAtmosphere_c``
* ``utils_c``

Development Environment Setup
------------------------------

For developers who want to contribute, additional tools are required:

.. code-block:: bash

   pip install sphinx sphinx_rtd_theme

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd docs
   make html

The built documentation can be viewed at ``docs/_build/html/index.html``.

Verification
------------

To verify that installation completed correctly, run the example:

.. code-block:: bash

   cd example
   python ../Trajectory_Optimization.py example-settings.json

If executed successfully, optimization results will be output.
