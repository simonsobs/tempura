#!/bin/bash

fflags="-qopenmp -fPIC"
f2pycomp="intelem"

# f90 shared library
target=lib${pylib}

# set directory
flibloc="../F90/"

pylibloc="../py/"
wrapdir="../wrap/"

