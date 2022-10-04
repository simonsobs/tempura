#!/bin/bash

# set directory
dir_scripts="scripts/"
dir_f90src="fortran/"
dir_pytempura="pytempura/"

# files to be scanned
scan="norm_lens.f90 norm_src.f90 norm_general.f90 noise_spec.f90"

# create python interface
python ${dir_scripts}/interface.py -moddir ${dir_f90src} -modname ${scan}
mv norm_lens.py norm_src.py norm_general.py noise_spec.py ${dir_pytempura}

