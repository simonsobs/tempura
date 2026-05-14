=======
tempura
=======

.. image:: https://readthedocs.org/projects/tempura/badge/?version=latest
           :target: https://tempura.readthedocs.io/en/latest/?badge=latest
		   :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/simonsobs/tempura/badge.svg?branch=master
		   :target: https://coveralls.io/github/simonsobs/tempura?branch=master

.. image:: https://badge.fury.io/py/pytempura.svg
		       :target: https://badge.fury.io/py/pytempura


	   
(T)ool for (E)fficient co(MPU)tation of mode-coupling estimato(R) norm(A)lization
---------------------------------------------------------------------------------

.. image:: tests/data/image.png

This package contains a python module to compute analytic normalization of quadratic estimators for lensing, cosmic birefringence, patchy tau and point sources, based on separable formula. 

The code was verified in the following studies:

* **Lensing Reconstruction, Delensing** \
  Namikawa & Nagata JCAP 09 (2014) 009, https://arxiv.org/abs/1405.6568
* **Cosmic Birefringence** \
  Namikawa et al. PRD 101 (2020) 083527, https://arxiv.org/abs/2001.10465
* **Patchy Reionization** \
  Namikawa PRD, 97 (2018) 063505, https://arxiv.org/abs/1711.00058

Installing
----------

Make sure your ``pip`` tool is up-to-date. To install ``pytempura``, run:

.. code-block:: console
		
   $ pip install pytempura --user

This will install a pre-compiled binary suitable for your system (only Linux and Mac OS X with Python>=3.10 are supported). 

If you require more control over your installation, e.g. using Intel compilers, please see the section below on compiling from source.

  
Compiling from source (advanced / development workflow)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to install from source is to use the ``pip`` tool,
with the ``--no-binary`` flag. This will download the source distribution
and compile it for you. Don't forget to make sure you have CC and FC set
if you have any problems.

For all other cases, below are general instructions.

First, download the source distribution or ``git clone`` this repository. You
can work from ``master`` or checkout one of the released version tags (see the
Releases section on Github). Then change into the cloned/source directory.

Once downloaded, you can install using ``pip install .`` inside the project
directory. We use the ``meson`` build system, which should be understood by
``pip`` (it will build in an isolated environment).

We suggest you then test the installation by running the unit tests. You
can do this by running ``pytest``.

To run an editable install, you will need to do so in a way that does not
have build isolation (as the backend build system, `meson` and `ninja`, actually
perform micro-builds on usage in this case):

.. code-block:: console
   
   $ pip install --upgrade pip meson ninja meson-python cython numpy
   $ pip install  --no-build-isolation --editable .

   
Examples
--------

You can find example codes at "tests" directory. 


Contact
-------

* Toshiya Namikawa (namikawa.toshiya9@gmail.com)
* Mathew Madhavacheril (mathewsyriac@gmail.com)


