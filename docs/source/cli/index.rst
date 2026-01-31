Command Line Interface
======================

Gelex provides a comprehensive command-line interface for genomic analysis. The tool follows a subcommand-based structure similar to ``git`` or ``plink2``.

Available Subcommands
---------------------

The following subcommands cover the complete genomic analysis workflow:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Subcommand
     - Description
   * - :doc:`fit`
     - Fit Bayesian models (BayesAlphabet) and estimate marker effects.
   * - :doc:`assoc`
     - Perform GWAS using mixed linear models (GBLUP) with LOCO.
   * - :doc:`grm`
     - Compute Genomic Relationship Matrices (Yang, Zeng, Vitezica).
   * - :doc:`predict`
     - Predict phenotypes for new samples using trained effects.
   * - :doc:`simulate`
     - Simulate phenotypes based on real genotypes and genetic architectures.

.. toctree::
   :hidden:

   fit
   assoc
   grm
   predict
   simulate

Global Options
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Option
     - Description
   * - ``-h, --help``
     - Show help message and exit.
   * - ``-v, --version``
     - Print version information and exit.

Basic Usage Pattern
-------------------

To get help for a specific subcommand, use:

.. code-block:: bash

   gelex <subcommand> --help
