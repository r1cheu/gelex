Command Line Interface
======================

The ``gelex`` binary follows a subcommand-based structure similar to ``git`` or ``plink2``.

Available Subcommands
---------------------

The subcommands are grouped by the workflow they belong to.

Genomic Selection
~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Subcommand
     - Description
   * - :doc:`mcmc`
     - Fit Bayesian genomic-prediction models (BayesAlphabet) via MCMC and estimate marker effects.
   * - :doc:`post`
     - Summarize posterior diagnostics from MCMC sample files.
   * - :doc:`predict`
     - Predict phenotypes for new samples using fitted SNP effects.

Association & Variance Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Subcommand
     - Description
   * - :doc:`grm`
     - Build additive/dominance genomic relationship matrices from PLINK data. Prerequisite for ``assoc`` and ``reml``.
   * - :doc:`reml`
     - Estimate variance components with AI-REML.
   * - :doc:`assoc`
     - Run mixed-model association testing with optional LOCO.

.. toctree::
   :hidden:

   mcmc
   post
   predict
   grm
   reml
   assoc

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
