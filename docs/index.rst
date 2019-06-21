.. CathPy documentation master file, created by
   sphinx-quickstart on Fri Aug 24 16:46:39 2018.

Welcome to CathPy's documentation!
==================================

.. image:: https://readthedocs.org/projects/cathpy/badge/?version=latest
  :target: https://cathpy.readthedocs.io/en/latest/?badge=latest

.. image:: https://travis-ci.com/UCL/cathpy.svg?branch=master
  :target: https://travis-ci.com/UCL/cathpy

.. image:: https://codecov.io/gh/UCL/cathpy/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/UCL/cathpy


**CathPy** is a Bioinformatics library and toolkit written for `Python3 <https://python.org>`_.
It is written and maintained by the maintainers of `CATH <http://www.cathdb.info>`_, the Protein 
Structure Classification Database at `UCL <https://www.ucl.ac.uk/orengo-group>`_, UK.

The code was mainly written for use within the CATH group and covers a range of tasks from manipulating
protein 3D structures and sequence alignments. This software has been
released to the public domain since some of the functionality may also be useful to groups outside 
of the CATH group.

Getting Started
---------------

This software is released regularly on `PyPI <https://www.pypi.org>`_ with 
`semantic versioning <https://semver.org>`_. The latest stable 
release can be installed in your local Python environment with:

::

    pip install cathpy

You can update an existing installation with the latest version by adding the ``--update`` flag.

::

    pip install --update cathpy

Note: in very old versions of ``pip`` you may find it easier to uninstall (``pip uninstall cathpy``) then reinstall 
(``pip install cathpy``).

Guide
-----

.. toctree::
    :maxdepth: 2

    tools
    api
    help

Code
----

The code is maintained on GitHub:

::

    git clone https://github.com/UCL/cathpy

Please submit any feature requests and bug reports as `GitHub Issues <https://github.com/UCL/cathpy/issues>`_ 
(accompanying Pull Requests are welcome!).


Contributors
------------

- Ian Sillitoe `@sillitoe <https://github.com/sillitoe>`_
- Sayoni Das `@sayonidas03 <https://github.com/sayonidas03>`_

If you want to see your name here, then submit a pull request to the code base...


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
