project('bezier', 'c', 'fortran', version: '2024.6.21', meson_version: '>=0.60.0')

c_args = []
numpy_dep = dependency('numpy', required: true)
py_dep = find_installation('python3')
py_dep_version = py_dep.version()
# NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
numpy_nodepr_define = 'NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'
add_project_arguments('-D' + numpy_nodepr_define, language: 'c')


# Fortran sources
# Listing types.f90 first as other files depend on it.
fortran_sources = files(
  'src/fortran/types.f90',
  'src/fortran/status.f90',
  'src/fortran/helpers.f90',
  'src/fortran/quadpack.f90',
  'src/fortran/curve.f90',
  'src/fortran/curve_intersection.f90',
  'src/fortran/triangle.f90',
  'src/fortran/triangle_intersection.f90'
)

# Fortran library
libbezier_static_lib = static_library(
  'bezier',
  fortran_sources,
  include_directories: 'src/fortran/include',
  install: false, # Not installed system-wide, just used for Cython extension
)

# Cython extension
cython_sources = custom_target(
  '_speedup_cython_compilation',
  input: 'src/python/bezier/_speedup.pyx',
  output: '_speedup.c',
  command: [
    find_program('cython'), '--fast-fail', '-o', '@OUTPUT@', '@INPUT@',
    '--include-path', meson.current_source_dir() / 'src/python/bezier',
  ],
  depends: libbezier_static_lib, # Ensure libbezier is built first
)

pymod = python_installation()
_speedup_module = pymod.extension_module(
  '_speedup',
  cython_sources,
  c_args: c_args,
  include_directories: [numpy_dep.include_directories(), 'src/fortran/include'],
  link_with: [libbezier_static_lib],
  install: true,
  subdir: 'bezier',
)

# Install the Python module
python_installation = meson.get_compiler('python').find_installation()
python_installation.install_sources(
    files('src/python/bezier/__init__.py', 'src/python/bezier/_algebraic_intersection.py', 'src/python/bezier/_geometric_intersection.py', 'src/python/bezier/_intersection_helpers.py', 'src/python/bezier/_plot_helpers.py', 'src/python/bezier/curve.py', 'src/python/bezier/surface.py', 'src/python/bezier/triangle.py'),
    subdir: 'bezier'
)
