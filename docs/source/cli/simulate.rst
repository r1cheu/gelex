.. _simulate-command:

simulate
========

Simulate phenotypes based on genetic data and specified parameters.

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex simulate --bfile <bfile> --causal <causal_file> [OPTIONS]

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
   * - ``--causal``
     - (required)
     - Causal variants file (TSV: ID, optional effect)
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
   * - ``--seed``
     - -1
     - Random seed for reproducibility (-1 for time-based)

Examples
--------

.. code-block:: bash
   :caption: Basic Phenotype Simulation

   gelex simulate \
     -b genotypes \
     --causal causal_variants.txt \
     -o simulated_phenotypes

.. code-block:: bash
   :caption: Custom Heritability

   gelex simulate \
     -b genotypes \
     --causal causal_variants.txt \
     --h2 0.3 \
     -o simulated_h2_0.3

.. code-block:: bash
   :caption: Reproducible Simulation

   gelex simulate \
     -b genotypes \
     --causal causal_variants.txt \
     --seed 12345 \
     -o simulated_reproducible
