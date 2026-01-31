Installation
============

Gelex can be installed via package managers (recommended for users) or built from source (recommended for developers).

.. note::
   Currently, Gelex only supports **Linux x86_64** architectures. Support for other platforms is planned for future releases.

Binary Installation (Recommended)
---------------------------------

The easiest way to install Gelex is through `pixi <https://pixi.sh>`_ or `conda`.

Using Pixi
~~~~~~~~~~

Install Gelex globally using pixi:

.. code-block:: bash

   pixi g install -c conda-forge -c https://prefix.dev/gelex gelex

Using Conda
~~~~~~~~~~~

You can install Gelex from the ``prefix.dev`` channel:

.. code-block:: bash

   conda install -c conda-forge -c https://prefix.dev/gelex gelex

Build from Source
-----------------

If you want to contribute to development or use the latest features, you can build Gelex from source.

Prerequisites
~~~~~~~~~~~~~

- **Pixi**: We use Pixi to manage dependencies and build environments. Install it from `pixi.sh <https://pixi.sh>`_.
- **C++ Compiler**: A compiler with C++23 support (e.g., GCC 13+, Clang 16+).

Build Steps
~~~~~~~~~~~

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/r1cheu/gelex.git
      cd gelex

2. Install dependencies and build:

   .. code-block:: bash

      # Install all dependencies via pixi
      pixi install

      # Build the debug version (includes tests)
      pixi run build-debug

      # Build the release version (optimized)
      pixi run build-release

3. Install to your local system:

   .. code-block:: bash

      # Install the release binary to ~/.local/bin
      pixi run install-release

Verification
------------

After installation, verify that Gelex is working correctly:

.. code-block:: bash

   gelex --help

