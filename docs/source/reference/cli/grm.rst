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

   gelex grm --bfile <genotype_prefix> [--mode <A|D|AD>] [--loco] [OPTIONS]

Required input is PLINK genotype prefix (``--bfile``).

Method Selection
----------------

Choose ``--geno-method`` based on scaling strategy.

Detailed formulas and method definitions: :ref:`genotype-processor-methods`.

.. list-table::
   :header-rows: 1
   :widths: 24 43 33

   * - Code (long name)
     - Use when
     - Notes
   * - ``OSH`` (``OrthStandardizeHWE``), default
     - You want the default orthogonal HWE-standardized GRM.
     - Orthogonal dominance, HWE moments. Best for most workflows.
   * - ``SH`` (``StandardizeHWE``)
     - You want HWE standardization without orthogonal dominance.
     - Simpler encoding, HWE moments.
   * - ``CH`` (``CenterHWE``)
     - You prefer HWE centering without variance scaling.
     - Preserves original scale, HWE moments.
   * - ``OCH`` (``OrthCenterHWE``)
     - You need orthogonal HWE centering (no scaling).
     - Matches assoc default; HWE moments.
   * - ``S``, ``C``, ``OS``, ``OC``, ``NS``, ``NC``
     - You want sample-based statistics instead of HWE-based.
     - More data-dependent estimates.

If unsure, use the default (``OSH``).

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-o, --out`` ``grm``
   Output prefix for GRM files.

``--geno-method`` ``OSH``
   GRM coding method. Accepts codes ``SH``, ``CH``, ``OSH``, ``OCH``,
   ``S``, ``C``, ``OS``, ``OC``, ``NS``, ``NC``. See
   :ref:`genotype-processor-methods` for what each code means.

.. rubric:: Matrix Selection

``--mode`` ``A``
   Effect mode(s) to compute: ``A`` (additive), ``D`` (dominance), or
   ``AD`` (both additive and dominance).

``--loco`` ``false``
   Compute chromosome-wise LOCO GRMs.

.. rubric:: Performance Options

``-c, --chunk-size`` ``10000``
   Number of SNPs per chunk. Lower values reduce memory usage.

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads. Must be non-negative; ``0`` leaves the thread
   count to the runtime default.

Output Files
------------

Output naming depends on whether you request one or multiple matrices.

.. list-table::
   :header-rows: 1
   :widths: 34 28 38

   * - Scenario
     - File pattern
     - Notes
   * - Single mode (``--mode A`` or ``--mode D``), no LOCO
     - ``<out>.A.bin/.id`` or ``<out>.D.bin/.id``
     - The effect suffix (``.A``/``.D``) is always written.
   * - ``--mode AD``, no LOCO
     - ``<out>.A.bin/.id`` and ``<out>.D.bin/.id``
     - One file pair per effect mode.
   * - LOCO enabled
     - ``<out>.<add|dom>.chrNN.bin/.id``
     - One file pair per chromosome (and effect mode); ``NN`` is the
       zero-padded chromosome number (``chr01``, ``chr02``, ...).

File structure follows :ref:`grm-format`.

Warnings and Notes
------------------

.. note::

   If ``--mode`` is not set, Gelex defaults to the additive GRM (``--mode A``).

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
      --mode D \
      --geno-method OCH \
      -o my_grm_dom

.. code-block:: bash
   :caption: Additive and Dominance Together

   gelex grm \
      -b genotypes \
      --mode AD \
      --geno-method OSH \
      -o my_grm_both

.. code-block:: bash
   :caption: LOCO GRM for GWAS

   gelex grm \
      -b genotypes \
      --mode A \
      --loco \
      --geno-method S \
      -o my_grm_loco

See Also
--------

- :ref:`assoc-command` for GWAS using GRM inputs.
- :ref:`grm-format` for GRM binary format details.
