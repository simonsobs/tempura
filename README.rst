=======
tempura
=======

.. image:: https://travis-ci.org/simonsobs/tempura.svg?branch=master
           :target: https://travis-ci.org/simonsobs/tempura

.. image:: https://readthedocs.org/projects/tempura/badge/?version=latest
           :target: https://tempura.readthedocs.io/en/latest/?badge=latest
		   :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/simonsobs/tempura/badge.svg?branch=master
		   :target: https://coveralls.io/github/simonsobs/tempura?branch=master

.. image:: https://badge.fury.io/py/tempura.svg
		       :target: https://badge.fury.io/py/tempura


	   
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


Installation
------------

Run 

.. code-block:: console
		
   $ python setup.py build_ext -i

  
You can then 'import tempura' if this repository is in your PYTHONPATH.

We are working on making this package available on PyPI (for a pip install).


Examples
--------

You can find example codes at "tests" directory. 


Contact
-------

* Toshiya Namikawa (namikawa.toshiya9@gmail.com)


