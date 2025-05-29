# Installing `bezier`

This guide provides instructions for installing the `bezier` Python package,
covering build dependencies, installation from PyPI or local source, and
how to run tests.

## 1. Build Dependencies

To build `bezier` from source, you'll need the following:

*   **Python:** Version >=3.10 (as specified in `pyproject.toml`).
*   **C Compiler:** A standard C compiler (e.g., GCC, Clang, MSVC).
*   **Fortran Compiler:** A Fortran compiler (e.g., gfortran).
*   **Build Tools:**
    *   `meson` - The Meson build system
    *   `ninja` - The Ninja build tool
    *   `cython` - For building Cython extensions
    *   `numpy` - Required for building the package
*   **Package Manager (optional but recommended):**
    *   `uv` - A modern Python package installer and resolver

**Installation of Dependencies:**

### Using uv (Recommended)

```bash
# Install uv if you haven't already
curl -sSf https://astral.sh/uv/install.sh | sh

# Create and activate a virtual environment
uv venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install build dependencies
uv pip install -U pip setuptools wheel
uv pip install meson-python cython numpy ninja
```

### Manual Installation

1. **Python:** Download from [python.org](https://www.python.org/) or use a version manager like `pyenv`.

2. **C/Fortran Compilers:**
   * **Linux (Debian/Ubuntu):** `sudo apt-get install gcc gfortran`
   * **Linux (Fedora):** `sudo dnf install gcc gfortran`
   * **macOS:** Install Xcode Command Line Tools and Homebrew, then: `brew install gcc`
   * **Windows:** Install MinGW-w64 for GCC/gfortran or use MSVC with Meson

3. **Build Tools:**
   ```bash
   pip install -U pip setuptools wheel
   pip install meson-python cython numpy ninja
   ```

## 2. Building the Package

The `bezier` package uses `meson-python` as its build backend, which integrates `meson` to compile the underlying Fortran library and Cython extensions.

### From PyPI (once published)

If a pre-compiled wheel is available for your platform and Python version, `pip` will install it directly:

```bash
pip install bezier
```

If a wheel is not available, `pip` will attempt to build the package from the source distribution. This will require the build dependencies listed above to be installed on your system.

### From Local Source (for development)

To install `bezier` from a local source tree (e.g., for development):

1. **Clone the repository and navigate to it:**
   ```bash
   git clone https://github.com/dhermes/bezier.git
   cd bezier
   ```

2. **Install in development mode (recommended):**
   This creates an "editable" installation where changes to the source code are immediately available without reinstalling.
   ```bash
   # Using uv (recommended)
   uv pip install --no-build-isolation -e .
   
   # Or using pip
   # pip install --no-build-isolation -e .
   ```

3. **For production installation (not recommended for development):**
   ```bash
   # Using uv
   uv pip install --no-build-isolation .
   
   # Or using pip
   # pip install --no-build-isolation .
   ```

   The `--no-build-isolation` flag ensures that the build uses the Python environment where you're installing the package, which can help avoid version conflicts.

## 3. Running Tests

Tests for `bezier` are managed using `nox`.

1.  **Install `nox`:**
    ```bash
    pip install nox
    ```

2.  **List available Nox sessions:**
    This command shows all defined test and utility sessions.
    ```bash
    nox -l
    ```

3.  **Run specific sessions:**
    You can run individual test suites or checks:
    ```bash
    # Run unit tests (uses default Python interpreter)
    nox -s unit

    # Run functional tests
    nox -s functional

    # Run doctests
    nox -s doctest

    # Run linting and code style checks
    nox -s lint
    ```
    You can also specify Python versions for sessions that support them, e.g., `nox -s "unit(python='3.10')"`.

4.  **Run all default sessions:**
    This typically includes main test suites and linting.
    ```bash
    nox
    ```

## 4. Build Troubleshooting

### Common Issues

1. **Missing Build Dependencies**
   - Ensure all build dependencies are installed (C/Fortran compilers, Meson, Ninja, Cython)
   - On Linux, you might need development packages like `python3-dev`

2. **Editable Installation Issues**
   - If you encounter issues with editable installs, try:
     ```bash
     pip uninstall -y bezier
     pip cache remove bezier
     pip install --no-build-isolation -e .
     ```

3. **Numpy Compatibility**
   - The package requires NumPy headers for building
   - If you get NumPy-related build errors, try:
     ```bash
     pip install -U numpy
     ```

## 5. Notes

*   **Build Process:** The package uses `meson-python` to build the Fortran library (`libbezier`) and Cython extensions as part of the standard Python package installation process.
*   **Build Configuration:** The build is configured by:
    - `pyproject.toml`: Specifies `meson-python` as the build backend and lists build dependencies
    - `meson.build`: Contains the Meson build script for the Fortran library and Cython extensions
    - `src/python/bezier/`: Contains the Python package source code
    - `src/fortran/`: Contains the Fortran source code

*   **Development Workflow:**
    - For development, use `pip install --no-build-isolation -e .`
    - The build system will automatically recompile changed Fortran/Cython files when the package is imported
    - To force a clean rebuild, delete the `build/` directory and reinstall
```
