.. _grm-command:

grm
===

Compute genomic relationship matrix (GRM) from PLINK binary files.

Use this command before :ref:`assoc-command` when you need mixed-model GWAS.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex grm -b genotypes -o my_grm

.. code-block:: bash
   :caption: Full Syntax Template

   gelex grm --bfile <genotype_prefix> [--add] [--dom] [--loco] [OPTIONS]

Required input is PLINK genotype prefix (``--bfile``).

Method Selection
----------------

Choose ``--geno-method`` based on scaling strategy.

Detailed formulas and method definitions: :ref:`genotype-processor-methods`.

.. list-table::
   :header-rows: 1
   :widths: 10 43 47

   * - Int
     - Use when
     - Notes
   * - ``3`` (default)
     - You want the default orthogonal HWE-standardized GRM.
     - Orthogonal dominance, HWE moments. Best for most workflows.
   * - ``1``
     - You want HWE standardization without orthogonal dominance.
     - Simpler encoding, HWE moments.
   * - ``2``
     - You prefer HWE centering without variance scaling.
     - Preserves original scale, HWE moments.
   * - ``4``
     - You need orthogonal HWE centering (no scaling).
     - Matches assoc default; HWE moments.
   * - ``5``–``8``
     - You want sample-based statistics instead of HWE-based.
     - More data-dependent estimates.

If unsure, use the default (``3``).

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-o, --out`` ``grm``
   Output prefix for GRM files.

``--geno-method`` ``3``
   GRM method as integer 1–8.
   ``1``\ =standardize-hwe, ``2``\ =center-hwe,
   ``3``\ =orth-standardize-hwe, ``4``\ =orth-center-hwe,
   ``5``\ =standardize, ``6``\ =center, ``7``\ =orth-standardize,
   ``8``\ =orth-center.
   See :ref:`genotype-processor-methods`.

.. rubric:: Matrix Selection

``--add`` ``false``
   Compute additive GRM.

``--dom`` ``false``
   Compute dominance GRM.

``--loco`` ``false``
   Compute chromosome-wise LOCO GRMs.

.. rubric:: Performance Options

``-c, --chunk-size`` ``10000``
   Number of SNPs per chunk. Lower values reduce memory usage.

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads (use ``-1`` for all cores).

Output Files
------------

Output naming depends on whether you request one or multiple matrices.

.. list-table::
   :header-rows: 1
   :widths: 34 28 38

   * - Scenario
     - File pattern
     - Notes
   * - Single matrix (additive or dominance), no LOCO
     - ``<out>.bin`` and ``<out>.id``
     - No ``.add``/``.dom`` suffix in this mode.
   * - Additive + dominance, no LOCO
     - ``<out>.add.bin/.id`` and ``<out>.dom.bin/.id``
     - One file pair per matrix type.
   * - LOCO enabled
     - ``<out>.<add|dom>.chrN.bin/.id``
     - One file pair per chromosome (and matrix type).

File structure follows :ref:`grm-format`.

Warnings and Notes
------------------

.. note::

   If neither ``--add`` nor ``--dom`` is set, Gelex defaults to additive GRM.

.. warning::

   ``--loco`` can generate many files for large chromosome sets. Ensure your
   downstream ``assoc --grm`` inputs use the matching LOCO prefix.

Examples
--------

.. code-block:: bash
   :caption: Default Additive GRM

   gelex grm \
      -b genotypes \
      -o my_grm

.. code-block:: bash
   :caption: Dominance GRM with Orthogonal HWE Centering

   gelex grm \
      -b genotypes \
      --dom \
      --geno-method 4 \
      -o my_grm_dom

.. code-block:: bash
   :caption: Additive and Dominance Together

   gelex grm \
      -b genotypes \
      --add \
      --dom \
      --geno-method 3 \
      -o my_grm_both

.. code-block:: bash
   :caption: LOCO GRM for GWAS

   gelex grm \
      -b genotypes \
      --add \
      --loco \
      --geno-method 5 \
      -o my_grm_loco

See Also
--------

- :ref:`assoc-command` for GWAS using GRM inputs.
- :ref:`grm-format` for GRM binary format details.
