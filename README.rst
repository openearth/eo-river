==============
EO River Tools
==============


.. image:: https://img.shields.io/pypi/v/eo-river.svg
        :target: https://pypi.python.org/pypi/eo-river

.. image:: https://img.shields.io/travis/openearth/eo-river.svg
        :target: https://travis-ci.org/openearth/eo-river

.. image:: https://readthedocs.org/projects/eo-river/badge/?version=latest
        :target: https://eo-river.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Earth Observation River Tools


* Free software: MIT license
* Documentation: https://eo-river.readthedocs.io.

Development
-----------

To install a new version, use the following (under Unix):

```bash
make bump
git push
make release
```

Then, manually release a new version on GitHub.

PyInstaller
-----------

To build exe on Windows, we use PyInstaller. Unfortunately, with the latest version of Anaconda it generates relatively large file (~280Mb) dut to MKL libraries.
The following workaround is used to avoid this: https://stackoverflow.com/questions/43886822/pyinstaller-with-pandas-creates-over-500-mb-exe/48846546#48846546.
This should not be too critical as we do not require a very high performance.

First, run the following (once):

```bash
conda create -n exe python=3
activate exe
pip install pandas pyinstaller pypiwin32 Click hydroengine geojson
echo hiddenimports = ['pandas._libs.tslibs.timedeltas'] > %CONDA_PREFIX%\Lib\site-packages\PyInstaller\hooks\hook-pandas.py
```
   
Then, the following command builds a new exe in EXE build\exe\dist\ directory:
```bash
scripts\build_exe.cmd
```


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
