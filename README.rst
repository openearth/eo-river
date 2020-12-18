========
EO River
========

A tool to extract water masks, network (skeleton) and water channel properties along network edges.

.. image:: https://user-images.githubusercontent.com/169821/47001405-a4e36e00-d12a-11e8-8de7-279f969e3f4e.png

.. _image:: https://img.shields.io/pypi/v/eo-river.svg
        :target: https://pypi.python.org/pypi/eo-river

.. _image:: https://img.shields.io/travis/openearth/eo-river.svg
        :target: https://travis-ci.org/openearth/eo-river

.. _image:: https://readthedocs.org/projects/eo-river/badge/?version=latest
        :target: https://eo-river.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status



Earth Observation River Tools


* Free software: MIT license
* Documentation: https://eo-river.readthedocs.io.

Development
-----------

To install a new version, use the following (under Unix):

  make bump
  
  git push
  
  make release


Then, manually release a new version on GitHub.

PyInstaller
-----------

To build exe on Windows, we use PyInstaller. Unfortunately, with the latest version of Anaconda it generates relatively large file (~280Mb) dut to MKL libraries. This should not be too critical as we do not require a very high performance.

First, run the following (one time only):

  conda create -n exe python=3
  
  activate exe
  
  pip install pandas pyinstaller pypiwin32 Click hydroengine geojson
  
  echo hiddenimports = ['pandas._libs.tslibs.timedeltas'] > %CONDA_PREFIX%/Lib/site-packages/PyInstaller/hooks/hook-pandas.py
   
Then, the following command builds a new exe in EXE build/exe/dist/ directory:

  scripts/build_exe.cmd


Features
--------

Usage:

  eo-river get-water-mask --region=extent.json --output=water_mask.json ... --smoothen

  eo-river generate-network --water-mask=water_mask.json --output=network.json ... --algorithm=mda|...

  eo-river generate-curvilinear-grid --network=network.json --output=grid.json ...

  eo-river generate-bathymetry --network=network.json --grid=grid.json --output=bath.nc ...

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
