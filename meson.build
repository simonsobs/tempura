# meson.build  (top level)

project(
  'pytempura',                       # project name
  ['c', 'fortran'],                  # languages we are going to compile
  default_options : [
    'warning_level=1',
    'buildtype=release'
  ],
)

py = import('python').find_installation(pure: false)
py_dep    = py.dependency()          # libpython + include path

# --- NumPy include dir -------------------------------------------------
numpy_inc = run_command(
  py, '-c', 'import numpy,sys;print(numpy.get_include())'
).stdout().strip()
inc_numpy = include_directories(numpy_inc)

# --- locate fortranobject.c correctly ---------------------------------
fortranobject_c_path = run_command(
  py, '-c',
  '''
import importlib.util, pathlib, sys, numpy
spec = importlib.util.find_spec("numpy.f2py")
print(pathlib.Path(spec.origin).parent / "src" / "fortranobject.c")
'''
).stdout().strip()

fortranobject_c = files(fortranobject_c_path)


f2py_src_inc = include_directories(
  run_command(
    py, '-c',
    '''
import importlib.util, pathlib
spec = importlib.util.find_spec("numpy.f2py")
print(pathlib.Path(spec.origin).parent / "src")
'''
  ).stdout().strip()
)


# ----------  common compile / link flags ----------
openmp_dep = dependency('openmp', required : true)

# Allow users to override optimisation flags like the old FCFLAGS
fc_env = run_command('sh', '-c', 'echo ${FCFLAGS:-}').stdout().strip().split()
if fc_env.length() == 0
  fc_env = ['-O3', '-fPIC']
endif

extra_fortran_args = [
  '-fopenmp', '-Wno-conversion', '-Wno-tabs', '-fPIC',
] + fc_env

add_project_arguments(extra_fortran_args, language : 'fortran')
add_project_link_arguments('-fopenmp', language : ['c','fortran'])

# ----------  libalkernel static library ----------
alkernel_src = files(
  'fortran/alkernel.f90',
  'fortran/norm_quad.f90',
)

alkernel_lib = static_library(
  'alkernel',
  alkernel_src,
  dependencies : openmp_dep,
  install      : false,          # it is only for internal linking
)

alkernel_dep = declare_dependency(link_with : alkernel_lib)

# ----------  _libtempura extension  ----------
libtempura_sources = files(
  'fortran/norm_lens.f90',
  'fortran/norm_src.f90',
  'fortran/norm_general.f90',
  'fortran/noise_spec.f90',
)

# Tell f2py to generate the wrapper code for all *.f90 files in one shot
f2py_wrappers = custom_target(
  '_libtempura_f2py',
  input  : libtempura_sources,
  output : ['_libtempura-f2pywrappers2.f90', '_libtempuramodule.c'],
  command: [py, '-m', 'numpy.f2py', '-m', '_libtempura', '--lower'] + libtempura_sources,
  depend_files : libtempura_sources,
)

py.extension_module(
  '_libtempura',
  # order is important: f2py wrappers, raw sources, fortranobject helper
  sources             : libtempura_sources + f2py_wrappers + [fortranobject_c],
  include_directories : [inc_numpy, f2py_src_inc],  
  dependencies        : [py_dep, openmp_dep, alkernel_dep],
  link_args           : ['-fopenmp'],
  subdir              : 'pytempura',      # installs into pytempura/_libtempura.so
  install             : true,
)

subdir('pytempura')

message('\nConfiguration summary:')
message('  NumPy include dir : ' + numpy_inc)
message('  Fortran flags     : ' + ' '.join(extra_fortran_args))
message('  OpenMP enabled    : ' + openmp_dep.found().to_string())
