### Tool for Efficient coMPUtation of mode-coupling estimatoR normAlization (TEMPURA)

This package contains a python module to compute analytic normalization of quadratic estimators for lensing, cosmic birefringence, patchy tau and point sources, based on separable formula. 

The code was verified in the following stiduies:

  - **Lensing Reconstruction, Delensing** \
   Namikawa & Nagata JCAP 09 (2014) 009, https://arxiv.org/abs/1405.6568
  - **Cosmic Birefringence** \
   Namikawa et al. PRD 101 (2020) 083527, https://arxiv.org/abs/2001.10465
  - **Patchy Reionization** \
   Namikawa PRD, 97 (2018) 063505, https://arxiv.org/abs/1711.00058


### Installation

  [1] Compile Fortran public software

  Go to "F90/src/" 

  [2] Create Fortran wrapper

  Go to the root directory and type one of the following command:

    - ./MAKEALL.sh all       (for generating all modules)
  
  You will find modules under "wrap/"


### Examples

You can find the example codes at "example" directory. 


### Contact

  tn334 at cam.ac.uk


