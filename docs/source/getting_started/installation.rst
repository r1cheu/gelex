Installation
============

Gelex can be installed via package managers (recommended for users) or built from source (recommended for developers).

.. note::
   Currently, Gelex only supports **Linux x86_64** architectures. Support for other platforms is planned for future releases.

Binary Installation (Recommended)
---------------------------------

Packages are published to the ``gelex`` channel on `prefix.dev <https://prefix.dev>`_
and are not mirrored on conda-forge or anaconda.org. Both commands below add that channel explicitly;
conda-forge is still required, as it supplies the runtime dependencies.

Using Pixi
~~~~~~~~~~

Install Gelex globally using `pixi <https://pixi.sh>`_:

.. code-block:: bash

   pixi global install -c conda-forge -c https://prefix.dev/gelex gelex

Using Conda
~~~~~~~~~~~

.. code-block:: bash

   conda install -c conda-forge -c https://prefix.dev/gelex gelex

Build from Source
-----------------

If you want to contribute to development or use the latest features, you can build Gelex from source.

Prerequisites
~~~~~~~~~~~~~

- **Pixi**: We use Pixi to manage dependencies and build environments. Install it from `pixi.sh <https://pixi.sh>`_.
- **C++ Compiler**: A compiler with C++23 support. Pixi provisions the toolchain (GCC 15) as part of the environment, so no separate system compiler is required.

Build Steps
~~~~~~~~~~~

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/r1cheu/gelex.git
      cd gelex

2. Install dependencies:

   .. code-block:: bash

      # Install all dependencies via pixi
      pixi install

3. Build the CLI. The ``build-cli`` task builds the core library together with
   the ``gelex`` command-line binary:

   .. code-block:: bash

      # Debug build (includes tests) -> build/cli-debug/apps/gelex
      pixi run build-cli

      # Optimized release build -> build/cli-release/apps/gelex
      pixi run build-cli release

4. Make the binary available on your ``PATH``. There is no dedicated install
   task; copy or symlink the built binary into a directory on your ``PATH``:

   .. code-block:: bash

      ln -s "$(pwd)/build/cli-release/apps/gelex" ~/.local/bin/gelex

Verification
------------

After installation, verify that Gelex is working correctly:

.. code-block:: bash

   gelex --help

