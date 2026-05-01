FIDASIM Python Libraries
========================

In addition to hosting the FIDASIM fortran source code, this repo contains python libraries used by FIDASIM to validate and generate the required input files to run the executable.

These python libraries can be found at:

$FIDASIM_DIR/lib/python/

which contains the following packages:
- efit/
- fidasim/

However, these are not packaged for quick installation by default.
For this reason, on 2026-04-22, we have introduced a pyproject.toml file to make them installable
as a standard python package to enable:

* Integration with IPS (Integrated Plasma Simulator) wrappers
* Clean reuse in standalone preprocessing scripts
* Removal of PYTHONPATH-based workflows

In other words, one can "pip install" this python package into working environment.
this makes the package accessible globally within the environment and thus importable

Intended Use Cases
==================

* IPS-FIDASIM wrapper development
* Preprocessing for FIDASIM runs
* Standalone validation and data preparation
* Integration with CQL3D workflows

Installation
=============

1- create a clean conda environment:
conda create -n MY_CONDA_ENV python=3.10
conda activate MY_CONDA_ENV

2- install the package inside the new conda env:
pip install -e /path/to/FIDASIM/lib/python

This installs the fidasim and efit Python modules in editable mode in the current conda environment

Usage
=====

Ater installation, while inside the new conda environment, you can access the python modules.
This can be done as follows:

inside a python script you can write:

from fidasim.preprocessing import prefida
import efit

Dependencies
============

This package depends on:

* numpy
* h5py
* scipy (< 1.11)
* scikit-image
* matplotlib

These are automatically installed via pyproject.toml

SciPy compatibility
===================
The EFIT utilities rely on legacy functions (cumtrapz) which are not available in newer versions of SciPy.

For this reason, scipy < 1.11 is enforced

Development workflow
====================
Editable install

pip install -e /path/to/FIDASIM/lib/python

This allows you to modify the source code and immediately use the changes without reinstalling

Testing Installation
====================
To test if the editable installation worked, in the terminal type the following (while inside the new conda environment):

python -c "from fidasim.preprocessing import prefida; import efit; print('OK')"
