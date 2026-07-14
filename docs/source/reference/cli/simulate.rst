.. _simulate-command:

simulate
========

Simulate phenotypes from real genotype data and user-defined architecture.

Use this command to create controlled datasets for method validation.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex simulate -b genotypes --h2 0.5 --add-var 0.01 --add-n 1000 -o sim_data

.. code-block:: bash
   :caption: Full Syntax Template

   gelex simulate --bfile <genotype_prefix> --h2 <heritability> [OPTIONS]

Required input is the genotype prefix (``--bfile``); you must also set at
least one of ``--h2`` or ``--d2``.

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-o, --out`` ``sim.phen``
   Output prefix/path root for simulation outputs.

``--h2``
   Additive heritability proportion in (0, 1). No default; set at least one
   of ``--h2`` or ``--d2``.

``--d2``
   Dominance heritability proportion in (0, 1). No default.

.. rubric:: Effect Architecture

``--add-var``
   Additive effect-class variances (one or more non-negative values).

``--add-n``
   SNP counts for additive effect classes; must match ``--add-var`` length.

``--dom-var``
   Dominance effect-class variances (one or more non-negative values).

``--dom-n``
   SNP counts for dominance effect classes; must match ``--dom-var`` length.

``--dom-pos-prob``
   Probability that dominance effects are positive, in (0, 1). Requires
   ``--d2``.

.. rubric:: Model

``--geno-method`` ``OS``
   Genotype coding method. Accepts codes ``SH``, ``CH``, ``OSH``, ``OCH``,
   ``S``, ``C``, ``OS``, ``OC``, ``NS``, ``NC``. See
   :ref:`genotype-processor-methods` for what each code means.

.. rubric:: Randomness

``--seed`` ``42``
   Random seed for reproducibility.

Output Files
------------

Simulation writes phenotype and causal-effect outputs using the ``--out`` root.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File pattern
     - Contents
     - Notes
   * - ``<out>.phen``
     - Simulated phenotype table (FID, IID, phenotype)
     - Main output for downstream ``fit`` or ``assoc``
   * - ``<out>.causal``
     - Causal SNP effects (id, additive[, dominance])
     - Ground truth for benchmarking

Warnings and Notes
------------------

.. warning::

   Keep ``h2 + d2 < 1`` to leave residual variance positive.

.. warning::

   ``--add-var`` and ``--add-n`` must have the same number of entries.

.. note::

   Dominance classes (``--dom-var`` and ``--dom-n``) are only used when
   ``--d2`` is greater than 0.

Examples
--------

.. code-block:: bash
   :caption: Basic Phenotype Simulation

   gelex simulate \
      -b genotypes \
      --h2 0.5 \
      --add-var 0.01 \
      --add-n 1000 \
      -o sim_basic

Expected outputs: ``sim_basic.phen``, ``sim_basic.causal``.

.. code-block:: bash
   :caption: Custom Heritability with Dominance

   gelex simulate \
      -b genotypes \
      --h2 0.3 \
      --add-var 0.01 \
      --add-n 1000 \
      --d2 0.1 \
      --dom-var 0.01 \
      --dom-n 500 \
      --seed 2026 \
      -o sim_dom

.. code-block:: bash
   :caption: Mixture Additive Effects (BayesR-style)

   gelex simulate \
      -b genotypes \
      --add-var 0 0.0001 0.001 0.01 \
      --add-n 900 50 30 20 \
      --h2 0.5 \
      --seed 42 \
      -o sim_mix

.. code-block:: bash
   :caption: Additive + Dominance Mixture

   gelex simulate \
      -b genotypes \
      --h2 0.4 \
      --d2 0.2 \
      --add-var 0 0.001 0.01 \
      --add-n 850 100 50 \
      --dom-var 0 0.001 \
      --dom-n 950 50 \
      --seed 42 \
      -o sim_arch

See Also
--------

- :ref:`mcmc-command` for training models on simulated phenotypes.
- :ref:`assoc-command` for GWAS benchmarking with known causal effects.
