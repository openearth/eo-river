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


Features
--------


Usage:

  pip install eoriver

  eoriver get-water-mask --region sues.geojson --filter-bounds filter.geojson --output sues-water.geojson
  
  get-water-mask
  get-water-network
  get-water-network-properties
  

  eoriver generate-network --water-mask=water_mask.json --output=network.json ... --algorithm=mda|...

  
Development
-----------

To install a new version, use the following (under Unix):

  make bump
  
  git push
  
  make release


Then, manually release a new version on GitHub.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
