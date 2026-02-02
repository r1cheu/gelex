.. _simulate-command:

simulate
========

Simulate phenotypes based on genetic data and specified parameters.

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex simulate --bfile <bfile> [OPTIONS]

Options
-------

.. rubric:: Data Files

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-b, --bfile``
     - (required)
     - PLINK binary file prefix (.bed/.bim/.fam)
   * - ``-o, --out``
     - "sim.phen"
     - Output file prefix for simulated phenotypes

.. rubric:: Simulation Parameters

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--h2``
     - 0.5
     - Narrow-sense heritability (0.0 - 1.0)
   * - ``--d2``
     - 0.0
     - Dominance variance proportion (0.0 - 1.0, h2+d2 < 1)
   * - ``--add-var``
     -
     - Variances for additive effect classes
   * - ``--add-prop``
     -
     - Proportions for additive effect classes (must match ``--add-var`` length, sum to 1)
   * - ``--dom-var``
     -
     - Variances for dominance effect classes
   * - ``--dom-prop``
     -
     - Proportions for dominance effect classes (must match ``--dom-var`` length, sum to 1)
   * - ``--seed``
     - -1
     - Random seed for reproducibility (-1 for time-based)

Examples
--------

.. code-block:: bash
   :caption: Basic Phenotype Simulation

   gelex simulate \
     -b genotypes \
     -o simulated_phenotypes

.. code-block:: bash
   :caption: Custom Heritability with Dominance

   gelex simulate \
     -b genotypes \
     --h2 0.3 --d2 0.1 \
     --seed 42

.. code-block:: bash
   :caption: Mixture Normal Effect Sizes (BayesR-style)

   gelex simulate \
     -b genotypes \
     --add-var 0 0.01 0.001 0.0001 \
     --add-prop 0.90 0.003 0.103 0.894 \
     --h2 0.5 --seed 42
