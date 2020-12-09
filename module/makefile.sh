#!/bin/bash

pylib=tempura

source ../sh/setup.sh

# compile link
modd="-I${flibloc}/mod"
libd="-L${flibloc}/lib"
link="-lutils"

# files to be compiled
scan="norm_lens.f90 norm_rot.f90 norm_tau.f90 norm_src.f90"

source ../sh/compile.sh ${1} ${2}

