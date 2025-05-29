project('bezier', 'c', 'fortran', version: '2024.6.21', meson_version: '>=0.60.0')

# Import Python module
py_mod = import('python')
py = py_mod.find_installation('python3')
py_dep = py.dependency()
py_version = py.language_version()

# Get numpy include path
numpy_inc = run_command(py, ['-c', 'import numpy as np; print(np.get_include())'], check: true).stdout().strip()
# Convert to relative path if it's under the source directory
if numpy_inc.startswith(meson.current_source_dir())
    numpy_inc = numpy_inc.split(meson.current_source_dir())[1].strip('/')
endif
numpy_inc_dir = include_directories(numpy_inc)
numpy_dep = dependency('numpy', required: true)

# C compiler arguments
c_args = []
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
    '-I', meson.current_source_dir() / 'src/python/bezier',
  ],
  depends: libbezier_static_lib, # Ensure libbezier is built first
)

# Build the Cython extension
_speedup_module = py.extension_module(
  '_speedup',
  cython_sources,
  c_args: c_args,
  include_directories: [numpy_inc_dir, include_directories('src/fortran/include')],
  link_with: [libbezier_static_lib],
  dependencies: [py_dep, numpy_dep],
  install: true,
  subdir: 'bezier',
)

# Install Python package
python_sources = [
    'src/python/bezier/__init__.py',
    'src/python/bezier/__config__.py',
    'src/python/bezier/_base.py',
    'src/python/bezier/_curve_helpers.py',
    'src/python/bezier/_geometric_intersection.py',
    'src/python/bezier/_helpers.py',
    'src/python/bezier/_intersection_helpers.py',
    'src/python/bezier/_plot_helpers.py',
    'src/python/bezier/_symbolic.py',
    'src/python/bezier/_triangle_helpers.py',
    'src/python/bezier/_triangle_intersection.py',
    'src/python/bezier/curve.py',
    'src/python/bezier/curved_polygon.py',
    'src/python/bezier/triangle.py',
    'src/python/bezier/hazmat/__init__.py',
    'src/python/bezier/hazmat/algebraic_intersection.py',
    'src/python/bezier/hazmat/clipping.py',
    'src/python/bezier/hazmat/curve_helpers.py',
    'src/python/bezier/hazmat/geometric_intersection.py',
    'src/python/bezier/hazmat/helpers.py',
    'src/python/bezier/hazmat/intersection_helpers.py',
    'src/python/bezier/hazmat/triangle_helpers.py',
    'src/python/bezier/hazmat/triangle_intersection.py'
]

# Install Python package files
py.install_sources(
    python_sources,
    subdir: 'bezier',
    pure: false
)

# Install the hazmat subpackage
py.install_sources(
    files(
        'src/python/bezier/hazmat/__init__.py',
        'src/python/bezier/hazmat/algebraic_intersection.py',
        'src/python/bezier/hazmat/clipping.py',
        'src/python/bezier/hazmat/curve_helpers.py',
        'src/python/bezier/hazmat/geometric_intersection.py',
        'src/python/bezier/hazmat/helpers.py',
        'src/python/bezier/hazmat/intersection_helpers.py',
        'src/python/bezier/hazmat/triangle_helpers.py',
        'src/python/bezier/hazmat/triangle_intersection.py'
    ),
    subdir: 'bezier/hazmat',
    pure: false
)
