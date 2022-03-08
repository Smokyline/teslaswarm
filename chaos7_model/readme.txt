
Overview
========

ChaosMagPy is a simple Python package for evaluating the
`CHAOS-7 <http://www.spacecenter.dk/files/magnetic-models/CHAOS-7/>`_ geomagnetic
field model. To quickly get started, download a complete working example
including the latest model under the "Forward code" section.

Documentation
-------------

The documentation of the current release is available on
`Read the Docs <https://chaosmagpy.readthedocs.io/en/stable/>`_.

|doi| |docs| |license|

.. |docs| image:: https://readthedocs.org/projects/chaosmagpy/badge/?version=stable
   :target: https://chaosmagpy.readthedocs.io/en/stable/?badge=stable
   :alt: Documentation Status

.. |license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: license.html

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3352398.svg
   :target: https://doi.org/10.5281/zenodo.3352398

References
----------

To reference ChaosMagPy in publications, please cite the package itself

https://doi.org/10.5281/zenodo.3352398

and all three of the following:

Finlay, C.C., Olsen, N., Kotsiaros, S., Gillet, N. and Toeffner-Clausen, L. (2016),
Recent geomagnetic secular variation from Swarm and ground observatories
as estimated in the CHAOS-6 geomagnetic field model Earth Planets Space,
Vol 68, 112. doi: 10.1186/s40623-016-0486-1

Olsen, N., Luehr, H., Finlay, C.C., Sabaka, T. J., Michaelis, I., Rauberg, J. and Toeffner-Clausen, L. (2014),
The CHAOS-4 geomagnetic field model, Geophys. J. Int., Vol 197, 815-827,
doi: 10.1093/gji/ggu033.

Olsen, N.,  Luehr, H.,  Sabaka, T.J.,  Mandea, M. ,Rother, M., Toeffner-Clausen, L. and Choi, S. (2006),
CHAOS — a model of Earth's magnetic field derived from CHAMP, Ørsted, and SAC-C magnetic satellite data,
Geophys. J. Int., vol. 166 67-75

License
=======

MIT License

Copyright (c) 2019 Clemens Kloss

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
============

ChaosMagPy relies on the following:

* python>=3.6
* numpy
* scipy
* pandas
* cython
* h5py
* cartopy>=0.17
* matplotlib>=3
* cdflib
* hdf5storage

Specific installation steps using the conda/pip package managers are as follows:

1. Install packages with conda:

   >>> conda install python numpy scipy pandas cython cartopy matplotlib h5py

2. Install packages with pip:

   >>> pip install cdflib hdf5storage

3. Finally install ChaosMagPy either with pip from PyPI:

   >>> pip install chaosmagpy

   or, if you have downloaded the `package files <https://pypi.org/project/chaosmagpy/#files>`_
   to the current working directory, with:

   >>> pip install chaosmagpy-x.x-py3-none-any.whl

   or, alternatively

   >>> pip install chaosmagpy-x.x.tar.gz

   replacing ``x.x`` with the correct version.

Contents
========

The directory contains the files/directories:

1. "chaosmagpy-x.x*.tar.gz": pip installable archive of the chaosmagpy package
   (version x.x*)

2. "chaos_examples.py": executable Python script containing several examples
   that can be run by changing the examples in line 16, save and run in the
   command line:

   >>> python chaos_examples.py

   example 1: Calculate CHAOS model field predictions from input coordinates
              and time and output simple data file
   example 2: Calculate and plot residuals between CHAOS model and
              Swarm A data (from L1b MAG cdf data file, example from May 2014).
   example 3: Calculate core field and its time derivatives for specified times
              and radii and plot maps
   example 4: Calculate static (i.e. small-scale crustal) magnetic field and
              plot maps (may take a few seconds)
   example 5: Calculate timeseries of the magnetic field at a ground
              observatory and plot
   example 6: Calculate external and associated induced fields described in SM
              and GSM reference systems and plot maps

3. "data/CHAOS-7.mat": mat-file containing CHAOS-7 model

4. "SW_OPER_MAGA_LR_1B_20180801T000000_20180801T235959_PT15S.cdf":
   cdf-file containing Swarm A magnetic field data from August 1, 2018.

5. directory called "html" containing the built documentation as
   html-files. Open "index.html" in your browser to access the main site.


Clemens Kloss (ancklo@space.dtu.dk)


Changelog
=========

Version 0.2.1
-------------
| **Date:** November 20, 2019
| **Release:** v0.2.1

Bugfixes
^^^^^^^^
* Corrected function :func:`chaosmagpy.coordinate_utils.zenith_angle` which was
  computing the solar zenith angle from ``phi`` defined as the hour angle and
  NOT the geographic longitude. The hour angle is measure positive towards West
  and negative towards East.

Version 0.2
-----------
| **Date:** October 3, 2019
| **Release:** v0.2

Features
^^^^^^^^
* Updated RC-index file to recent version (August 2019)
* Added option ``nmin`` to :func:`chaosmagpy.model_utils.synth_values`.
* Vectorized :func:`chaosmagpy.data_utils.mjd2000`,
  :func:`chaosmagpy.data_utils.mjd_to_dyear` and
  :func:`chaosmagpy.data_utils.dyear_to_mjd`.
* New function :func:`chaosmagpy.coordinate_utils.local_time` for a simple
  computation of the local time.
* New function :func:`chaosmagpy.coordinate_utils.zenith_angle` for computing
  the solar zenith angle.
* New function :func:`chaosmagpy.coordinate_utils.gg_to_geo` and
  :func:`chaosmagpy.coordinate_utils.geo_to_gg` for transforming geodetic and
  geocentric coordinates.
* Added keyword ``start_date`` to
  :func:`chaosmagpy.coordinate_utils.rotate_gauss_fft`
* Improved performance of :meth:`chaosmagpy.chaos.CHAOS.synth_coeffs_sm` and
  :meth:`chaosmagpy.chaos.CHAOS.synth_coeffs_gsm`.
* Automatically import :func:`chaosmagpy.model_utils.synth_values`.

Deprecations
^^^^^^^^^^^^
* Rewrote :func:`chaosmagpy.data_utils.load_matfile`: now traverses matfile
  and outputs dictionary.
* Removed ``breaks_euler`` and ``coeffs_euler`` from
  :class:`chaosmagpy.chaos.CHAOS` class
  attributes. Euler angles are now handled as :class:`chaosmagpy.chaos.Base`
  class instance.

Bugfixes
^^^^^^^^
* Fixed collocation matrix for unordered collocation sites. Endpoint now
  correctly taken into account.

Version 0.1
-------------
| **Date:** May 10, 2019
| **Release:** v0.1

Features
^^^^^^^^
* New CHAOS class method :meth:`chaosmagpy.chaos.CHAOS.synth_euler_angles` to
  compute euler angles for the satellites from the CHAOS model (used to rotate
  vectors from magnetometer frame to the satellite frame).
* Added CHAOS class methods :meth:`chaosmagpy.chaos.CHAOS.synth_values_tdep`,
  :meth:`chaosmagpy.chaos.CHAOS.synth_values_static`,
  :meth:`chaosmagpy.chaos.CHAOS.synth_values_gsm` and
  :meth:`chaosmagpy.chaos.CHAOS.synth_values_sm` for field value computation.
* RC index file now stored in HDF5 format.
* Filepaths and other parameters are now handled by a configuration dictionary
  called ``chaosmagpy.basicConfig``.
* Added extrapolation keyword to the BaseModel class
  :meth:`chaosmagpy.chaos.Base.synth_coeffs`, linear by default.
* :func:`chaosmagpy.data_utils.mjd2000` now also accepts datetime class
  instances.
* :func:`chaosmagpy.data_utils.load_RC_datfile` downloads latest RC-index file
  from the website if no file is given.

Bugfixes
^^^^^^^^
* Resolved issue in :func:`chaosmagpy.model_utils.degree_correlation`.
* Changed the date conversion to include hours and seconds not just the day
  when plotting the timeseries.

Version 0.1a3
-------------
| **Date:** February 19, 2019
| **Release:** v0.1a3

Features
^^^^^^^^
* New CHAOS class method :meth:`chaosmagpy.chaos.CHAOS.save_matfile` to output
  MATLAB compatible files of the CHAOS model (using the ``hdf5storage``
  package).
* Added ``epoch`` keyword to basevector input arguments of GSM, SM and MAG
  coordinate systems.

Bugfixes
^^^^^^^^
* Fixed problem of the setup configuration for ``pip`` which caused importing
  the package to fail although installation was indicated as successful.

Version 0.1a2
-------------
| **Date:** January 26, 2019
| **Release:** v0.1a2

Features
^^^^^^^^
* :func:`chaosmagpy.data_utils.mjd_to_dyear` and
  :func:`chaosmagpy.data_utils.dyear_to_mjd` convert time with microseconds
  precision to prevent round-off errors in seconds.
* Time conversion now uses built-in ``calendar`` module to identify leap year.

Bugfixes
^^^^^^^^
* Fixed wrong package requirement that caused the installation of
  ChaosMagPy v0.1a1 to fail with ``pip``. If installation of v0.1a1 is needed,
  use ``pip install --no-deps chaosmagpy==0.1a1`` to ignore faulty
  requirements.


Version 0.1a1
-------------
| **Date:** January 5, 2019
| **Release:** v0.1a1

Features
^^^^^^^^
* Package now supports Matplotlib v3 and Cartopy v0.17.
* Loading shc-file now converts decimal year to ``mjd2000`` taking leap years
  into account by default.
* Moved ``mjd2000`` from ``coordinate_utils`` to ``data_utils``.
* Added function to compute degree correlation.
* Added functions to compute and plot the power spectrum.
* Added flexibility to the function synth_values: now supports NumPy
  broadcasting rules.
* Fixed CHAOS class method synth_coeffs_sm default source parameter: now
  defaults to ``'external'``.

Deprecations
^^^^^^^^^^^^
* Optional argument ``source`` when saving shc-file has been renamed to
  ``model``.
* ``plot_external_map`` has been renamed to ``plot_maps_external``
* ``synth_sm_field`` has been renamed to ``synth_coeffs_sm``
* ``synth_gsm_field`` has been renamed to ``synth_coeffs_gsm``
* ``plot_static_map`` has been renamed to ``plot_maps_static``
* ``synth_static_field`` has been renamed to ``synth_coeffs_static``
* ``plot_tdep_maps`` has been renamed to ``plot_maps_tdep``
* ``synth_tdep_field`` has been renamed to ``synth_coeffs_tdep``


Version 0.1a0
-------------
| **Date:** October 13, 2018
| **Release:** v0.1a0

Initial release to the users for testing.
