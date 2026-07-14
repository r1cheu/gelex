.. _post-command:

post
====

Summarize posterior diagnostics from MCMC sample files.

Use this command after a Bayesian run (see :ref:`mcmc-command`) to compute
posterior means, medians, standard deviations, highest posterior density
intervals (HPDIs), effective sample sizes (ESS), and Gelman-Rubin ``R-hat``
diagnostics from one or more ``.draws`` files. Passing multiple prefixes
enables cross-chain ``R-hat`` convergence diagnostics.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex post --in mcmc_run

.. code-block:: bash
   :caption: Full Syntax Template

   gelex post --in <prefix...> [OPTIONS]

The only required input is one or more MCMC output prefixes (``--in``); each
prefix is read as ``<prefix>.draws``.

Options
-------

.. rubric:: Quick Start Options

``--in`` ``required``
   One or more MCMC output prefixes; reads ``<prefix>.draws``. Provide
   several prefixes (separate chains) to obtain ``R-hat`` diagnostics.

``-g, --gfile``
   Genotype binary prefix for per-SNP diagnostics (for example, GEBV variance
   and heritability). Required only when SNP-level diagnostics are wanted.

``-o, --out`` ``gelex_post``
   Output prefix for the run log.

.. rubric:: Model Configuration

``--hdpi`` ``0.94``
   Highest posterior density interval (HPDI) probability mass. Must be in the
   open interval ``(0, 1)``.

Output Files
------------

The posterior diagnostic table is printed to the console. In addition:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File pattern
     - Contents
   * - ``<out>.log``
     - Run log capturing the printed MCMC summary (parameter mean, median, SD,
       HPDI, ESS, and ``R-hat``).

Examples
--------

.. code-block:: bash
   :caption: Single Chain Summary

   gelex post \
      --in mcmc_run \
      -o post_run

.. code-block:: bash
   :caption: Multi-Chain R-hat Diagnostics

   gelex post \
      --in chain1 chain2 chain3 \
      -o post_multi

.. code-block:: bash
   :caption: With SNP-Level Diagnostics

   gelex post \
      --in mcmc_run \
      --gfile genotypes \
      --hdpi 0.9 \
      -o post_snp

See Also
--------

- :ref:`mcmc-command` for producing ``.draws`` files.
