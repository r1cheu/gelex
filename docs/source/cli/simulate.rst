.. _simulate-command:

simulate
========

Simulate phenotypes from real genotype data and user-defined architecture.

Use this command to create controlled datasets for method validation.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex simulate -b genotypes -o sim_data

.. code-block:: bash
   :caption: Full Syntax Template

   gelex simulate --bfile <genotype_prefix> [OPTIONS]

Required input is genotype prefix (``--bfile``).

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-o, --out`` ``sim.phen``
   Output prefix/path root for simulation outputs.

``--h2`` ``0.5``
   Additive heritability proportion.

``--d2`` ``0.0``
   Dominance heritability proportion.

.. rubric:: Effect Architecture

``--add-var`` ``0.01``
   Additive effect-class variances (one or more values).

``--add-prop`` ``1.0``
   Additive effect-class proportions; must match ``--add-var`` length and sum to 1.

``--dom-var`` ``0.01``
   Dominance effect-class variances.

``--dom-prop`` ``1.0``
   Dominance effect-class proportions; must match ``--dom-var`` length and sum to 1.

``--intercept`` ``0.0``
   Mean term added to simulated phenotypes.

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
     - Causal SNP effects and class assignments
     - Ground truth for benchmarking

Warnings and Notes
------------------

.. warning::

   Keep ``h2 + d2 < 1`` to leave residual variance positive.

.. warning::

   ``--add-var`` and ``--add-prop`` must have the same number of entries, and
   additive proportions must sum to 1.

.. note::

   Dominance classes (``--dom-var`` and ``--dom-prop``) are only used when
   ``--d2`` is greater than 0.

Examples
--------

.. code-block:: bash
   :caption: Basic Phenotype Simulation

   gelex simulate \
      -b genotypes \
      -o sim_basic

Expected outputs: ``sim_basic.phen``, ``sim_basic.causal``.

.. code-block:: bash
   :caption: Custom Heritability with Dominance

   gelex simulate \
      -b genotypes \
      --h2 0.3 \
      --d2 0.1 \
      --seed 2026 \
      -o sim_dom

.. code-block:: bash
   :caption: Mixture Additive Effects (BayesR-style)

   gelex simulate \
      -b genotypes \
      --add-var 0 0.0001 0.001 0.01 \
      --add-prop 0.90 0.05 0.03 0.02 \
      --h2 0.5 \
      --seed 42 \
      -o sim_mix

.. code-block:: bash
   :caption: Additive + Dominance Mixture with Intercept

   gelex simulate \
      -b genotypes \
      --h2 0.4 \
      --d2 0.2 \
      --add-var 0 0.001 0.01 \
      --add-prop 0.85 0.10 0.05 \
      --dom-var 0 0.001 \
      --dom-prop 0.95 0.05 \
      --intercept 1.5 \
      --seed 42 \
      -o sim_arch

See Also
--------

- :ref:`fit-command` for training models on simulated phenotypes.
- :ref:`assoc-command` for GWAS benchmarking with known causal effects.
