.. _grm-command:

grm
===

Compute genomic relationship matrix (GRM) from PLINK binary files.

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex grm --bfile <bfile> [OPTIONS]

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
     - "grm"
     - Output file prefix

.. rubric:: GRM Options

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-m, --method``
     - yang
     - Computation method: ``yang``, ``su``, ``zeng``, ``vitezica``
   * - ``-c, --chunk-size``
     - 10000
     - Chunk size for memory-efficient computation
   * - ``-t, --threads``
     - 12
     - Number of threads (-1 for all cores)
   * - ``--add``
     - false
     - Compute additive GRM
   * - ``--dom``
     - false
     - Compute dominance GRM
   * - ``--loco``
     - false
     - Compute separate GRM for each chromosome

Examples
--------

.. code-block:: bash
   :caption: Standard Additive GRM (Yang method)

   gelex grm -b genotypes --add -o my_grm

.. code-block:: bash
   :caption: Alternative Method (Su et al.)

   gelex grm -b genotypes --add -m su -o my_grm_su

.. code-block:: bash
   :caption: Dominance GRM

   gelex grm -b genotypes --dom -o my_grm_dom

.. code-block:: bash
   :caption: LOCO GRM Calculation

   gelex grm -b genotypes --add --loco -o my_grm_loco

.. code-block:: bash
   :caption: High Performance (Multi-threading)

   gelex grm -b genotypes --add --threads 16 -o my_grm_fast
