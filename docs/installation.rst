Installation
============

Requirements
------------

GELATO requires the following Python libraries:

* Python 3.8 or higher
* NumPy
* SciPy
* pandas
* pybind11
* simplekml
* pyoptsparse
* C++ compiler (for building PyBind11 modules)

**NLP Solver (Required)**

GELATO requires either IPOPT or SNOPT as a nonlinear programming (NLP) solver:

* **IPOPT** (recommended): Open-source interior point optimizer
* **SNOPT** (commercial): Sparse nonlinear optimizer

.. note::
   By using conda package manager, you can automatically install and set up IPOPT with pyoptsparse.

Installation Steps
------------------

Method 1: Using Conda (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install GELATO with all dependencies including IPOPT:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/istellartech/GELATO.git
   cd GELATO

   # Create conda environment and install dependencies with IPOPT
   conda install -c conda-forge numpy scipy pandas pybind11 simplekml
   conda install -c conda-forge ipopt pyoptsparse

   # Build C++ extension modules
   make

Method 2: Manual Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Clone the Repository
"""""""""""""""""""""""

.. code-block:: bash

   git clone https://github.com/istellartech/GELATO.git
   cd GELATO

2. Install Python Dependencies
"""""""""""""""""""""""""""""""

.. code-block:: bash

   pip install numpy scipy pandas pybind11 simplekml

3. Install NLP Solver
""""""""""""""""""""""

**Option A: Install IPOPT via conda**

.. code-block:: bash

   conda install -c conda-forge ipopt pyoptsparse

**Option B: Build pyoptsparse from source**

If you want to use a self-built IPOPT or SNOPT, you must build pyoptsparse from source.
Please refer to the `pyoptsparse installation guide <https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/optimizers/IPOPT.html>`_ for detailed instructions.

.. code-block:: bash

   # Install pyoptsparse from source
   git clone https://github.com/mdolab/pyoptsparse.git
   cd pyoptsparse
   pip install .

4. Build C++ Extension Modules
"""""""""""""""""""""""""""""""

GELATO uses C++ modules bound with PyBind11 for performance-critical calculations.
Build them with the following command:

.. code-block:: bash

   make

This will compile the following C++ extension modules:

* ``coordinate_c``: Coordinate transformation functions
* ``dynamics_c``: Rocket dynamics calculations
* ``IIP_c``: Instantaneous Impact Point calculations
* ``USStandardAtmosphere_c``: Atmospheric model
* ``utils_c``: Utility functions

Verification
------------

To verify that the installation completed correctly, run the example:

.. code-block:: bash

   cd example
   python ../Trajectory_Optimization.py example-settings.json

If executed successfully, optimization results will be output to a CSV file.

Troubleshooting
---------------

**Build Errors**

If you encounter errors during ``make``, ensure you have a C++ compiler installed:

* **Linux**: ``sudo apt-get install build-essential``
* **macOS**: Install Xcode Command Line Tools: ``xcode-select --install``
* **Windows**: Install Visual Studio with C++ support

**IPOPT/SNOPT Issues**

If pyoptsparse cannot find IPOPT or SNOPT:

1. Verify installation: ``python -c "from pyoptsparse import IPOPT"``
2. If using conda, ensure the environment is activated
3. Check the `pyoptsparse documentation <https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/>`_ for troubleshooting

Development Environment Setup
------------------------------

For developers who want to contribute to GELATO, additional tools are required:

.. code-block:: bash

   pip install sphinx sphinx_rtd_theme

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

To build the documentation locally:

.. code-block:: bash

   cd docs
   make html

The built HTML documentation can be viewed at ``docs/_build/html/index.html``.
