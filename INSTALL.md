# Installing `bezier`

This guide provides instructions for installing the `bezier` Python package,
covering build dependencies, installation from PyPI or local source, and
how to run tests.

## 1. Build Dependencies

To build `bezier` from source, you'll need the following:

*   **Python:** Version >=3.10 (as specified in `pyproject.toml`).
*   **C Compiler:** A standard C compiler (e.g., GCC, Clang, MSVC).
*   **Fortran Compiler:** A Fortran compiler (e.g., gfortran).
*   **Meson Build System:** The `meson` build system.
*   **Ninja:** The `ninja` build tool, which Meson typically uses.

**Installation of Dependencies:**

*   **Python:** Download from [python.org](https://www.python.org/) or use a version manager like `pyenv`.
*   **C/Fortran Compilers:**
    *   **Linux (Debian/Ubuntu):** `sudo apt-get install gcc gfortran`
    *   **Linux (Fedora):** `sudo dnf install gcc gfortran`
    *   **macOS:** Apple Command Line Tools (includes Clang) usually provide GCC and gfortran can be installed via Homebrew: `brew install gcc`.
    *   **Windows:** MinGW-w64 (for GCC/gfortran) is a common choice. MSVC can be used for C if Meson is configured accordingly.
*   **Meson and Ninja:**
    *   It's recommended to install these via `pip`: `pip install meson ninja`
    *   Alternatively, follow the official installation guides:
        *   [Meson Installation](https://mesonbuild.com/Getting-meson.html)
        *   [Ninja Installation](https://ninja-build.org/manual.html#_installation)

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

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/dhermes/bezier.git
    ```

2.  **Navigate to the directory:**
    ```bash
    cd bezier
    ```

3.  **Install using `pip`:**
    This command will build and install the package. `meson-python` handles the invocation of `meson` and `ninja` automatically.
    ```bash
    pip install .
    ```

4.  **Install in editable mode (for development):**
    This allows you to make changes to the source code and have them reflected immediately without reinstalling.
    ```bash
    pip install -e .
    ```

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

## 4. Notes

*   **Simplified Build Process:** With the transition to `meson-python`, the Fortran library (`libbezier`) and Cython extensions are compiled and linked as part of the standard Python package installation process (e.g., `pip install .`). This simplifies the previous multi-step build that required separate CMake configuration and installation of `libbezier`.
*   **Build Configuration:** The build process is now primarily defined by `pyproject.toml` (which specifies `meson-python` as the backend and lists build dependencies) and `meson.build` (which contains the Meson build script for the Fortran library and Cython extensions).
```
