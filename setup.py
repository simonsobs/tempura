#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
from __future__ import print_function
import setuptools
from setuptools import find_packages
from distutils.errors import DistutilsError
from numpy.distutils.core import setup, Extension, build_ext, build_src
import versioneer
import os, sys
import subprocess as sp
import numpy as np
import glob,shutil
build_ext = build_ext.build_ext
build_src = build_src.build_src


compile_opts = {
    'extra_f90_compile_args': ['-fopenmp', '-Wno-conversion', '-Wno-tabs', '-fPIC'],
    }

# Set compiler options
# Windows
if sys.platform == 'win32':
    raise DistUtilsError('Windows is not supported.')
# Mac OS X - needs gcc (usually via HomeBrew) because the default compiler LLVM (clang) does not support OpenMP
#          - with gcc -fopenmp option implies -pthread
elif sys.platform == 'darwin':
    try:
        sp.check_call('scripts/osx.sh', shell=True)
    except sp.CalledProcessError:
        raise DistutilsError('Failed to prepare Mac OS X properly. See earlier errors.')
    gccpath = glob.glob('/usr/local/bin/gcc-*')
    if gccpath:
        # Use newest gcc found
        sint = lambda x: int(x) if x.isdigit() else 0
        gversion = str(max([sint(os.path.basename(x).split('-')[1]) for x in gccpath]))
        os.environ['CC'] = 'gcc-' + gversion
        os.environ['CXX'] = os.environ['CC'].replace("gcc","g++")
        os.environ['FC'] = os.environ['CC'].replace("gcc","gfortran")
        rpath = '/usr/local/opt/gcc/lib/gcc/' + gversion + '/'
    else:
        os.system("which gcc")
        os.system("find / -name \'gcc\'")
        raise Exception('Cannot find gcc in /usr/local/bin. tempura requires gcc to be installed - easily done through the Homebrew package manager (http://brew.sh). Note: gcc with OpenMP support is required.')
    compile_opts['extra_link_args'] = ['-fopenmp', '-Wl,-rpath,' + rpath]
# Linux
elif sys.platform == 'linux':
    compile_opts['extra_link_args'] = ['-fopenmp']


def pip_install(package):
    import pip
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])

# with open('README.rst') as readme_file:
#     readme = readme_file.read()

# with open('HISTORY.rst') as history_file:
#     history = history_file.read()

requirements =  ['numpy>=1.16',
                 'setuptools>=39',
                 'matplotlib>=2.0',
                 'pytest-cov>=2.6',
                 'coveralls>=1.5',
                 'pytest>=4.6']


test_requirements = ['pip>=9.0',
                     'bumpversion>=0.5.',
                     'matplotlib>=2.0',
                     'wheel>=0.30',
                     'watchdog>=0.8',
                     'flake8>=3.5',
                     'coverage>=4.5',
                     'Sphinx>=1.7',
                     'twine>=1.10',
                     'numpy>=1.16',
                     'setuptools>=39.2',
                     'pytest-cov>=2.6',
                     'coveralls>=1.5',
                     'pytest>=4.6']
    
    
    

fcflags = os.getenv('FCFLAGS')
if fcflags is None or fcflags.strip() == '':
    fcflags = ['-O3','-fPIC']
else:
    print('User supplied fortran flags: ', fcflags)
    print('These will supersede other optimization flags.')
    fcflags = fcflags.split()
    
compile_opts['extra_f90_compile_args'].extend(fcflags)
compile_opts['extra_f77_compile_args'] = compile_opts['extra_f90_compile_args']

#snames = ['norm_lens']
#snames = ['norm_lens','norm_src','norm_rot','norm_tau']
snames = ['norm_lens','norm_src','norm_general']

class CustomBuild(build_ext):
    def run(self):
        print("Running build...")
        ret = build_ext.run(self)
        fs = glob.glob(f'_libtempura*.so')
        try:
            assert len(fs)==1
        except:
            print(fs)
            raise
        shutil.move(fs[0],f'pytempura/{fs[0]}')
        return ret

# Cascade your overrides here.
cmdclass = {
    'build_ext': CustomBuild,
}
cmdclass = versioneer.get_cmdclass(cmdclass)


setup(
    author="Toshiya Namikawa",
    author_email='namikawa.toshiya9@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Normalization for mode-coupling estimators",
    package_dir={"pytempura": "pytempura"},
    entry_points={
    },
    ext_modules=[
        Extension('_libtempura',
                  sources=[f'fortran/{sname}.f90' for sname in snames],
                  libraries = ['alkernel'],
                  #libraries = ['alkernel','norm_quad'],
                  **compile_opts),
    ],
    #libraries = [('alkernel', dict(sources=['fortran/alkernel.f90'],**compile_opts)),
    #    ('norm_quad', dict(sources=['fortran/norm_quad.f90'],**compile_opts))],
    #libraries = [('alkernel', dict(sources=['fortran/alkernel.f90'],**compile_opts))],
    libraries = [('alkernel', dict(sources=['fortran/alkernel.f90','fortran/norm_quad.f90'],**compile_opts))],
    install_requires=requirements,
    license="BSD license",
    # long_description=readme + '\n\n' + history,
    include_package_data=True,    
    keywords='CMB lensing',
    name='pytempura',
    packages=find_packages(),
    test_suite='pytempura.tests',
    tests_require=test_requirements,
    url='https://github.com/simonsobs/tempura',
    version=versioneer.get_version(),
    zip_safe=False,
    cmdclass=cmdclass,
    # scripts=['scripts/test-tempura']
)

print('\n[setup.py request was successful.]')

os.system('./scripts/generate_interface.sh')

print('\n[python interface files are created under pytempura/]')

