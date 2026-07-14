.. _assoc-command:

assoc
=====

Perform genome-wide association study (GWAS) using mixed linear models.

Use this command with GRM input from :ref:`grm-command`.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex assoc -b genotypes -p phenotypes.tsv --grm my_grm -o gwas_run

.. code-block:: bash
   :caption: Full Syntax Template

   gelex assoc --pheno <pheno_file> --bfile <genotype_prefix> --grm <grm_prefix...> [OPTIONS]

Required inputs are the phenotype file (``--pheno``), the genotype prefix
(``--bfile``), and at least one random effect. A random effect can be a GRM
(``--grm``), a discrete random effect (``--drand``), a quantitative random
effect (``--qrand``), or an interaction (``--interaction``). ``--grm`` is the
most common choice and is no longer mandatory.

Method Selection
----------------

Pick model and preprocessing strategy before tuning runtime options.

.. list-table::
   :header-rows: 1
   :widths: 22 43 35

   * - Option
     - Use when
     - Trade-off
   * - ``--mode A``
     - You are running a standard additive-effect GWAS (Wald test, df=1).
     - Fast and robust default for most analyses.
   * - ``--mode D``
     - You want to test dominance effects (single Wald test, df=1).
     - Requires compatible GRMs and often larger sample sizes.
   * - ``--mode AD``
     - You want a joint additive + dominance test (df=2).
     - Tests additive and dominance signals together.
   * - ``--transform none``
     - The phenotype is already approximately normal.
     - Keeps interpretation on the original trait scale.
   * - ``--transform dint`` / ``iint``
     - The phenotype distribution is skewed or heavy-tailed.
     - Often improves calibration, but effect sizes are on transformed scale.
   * - ``--geno-method OCH``
     - You want the default genotype preprocessing pipeline.
     - Orthogonal HWE centering. Good default for stability and comparability.

.. warning::

   If you use ``--transform`` with ``--mode D``, dominance signals may be
   attenuated, which can reduce power to detect dominance effects.


Options
-------

.. rubric:: Quick Start Options

``-p, --pheno`` ``required``
   Phenotype TSV file in format ``FID IID trait1 ...``.

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-o, --out`` ``gelex``
   Output prefix for GWAS results.

.. rubric:: Random Effect Inputs

At least one random effect is required; ``--grm`` is the usual choice.

``--grm``
   One or more GRM prefixes (``<prefix>.bin/.id``). Each prefix contributes one
   variance component.

``--drand``
   Discrete random-effect TSV (``FID IID factor1 ...``); each factor column
   becomes a variance component via a one-hot ``ZZ^T`` kernel.

``--qrand``
   One or more quantitative random-effect matrix TSVs (``FID IID value1 ...``);
   each file forms one linear-kernel ``ZZ^T`` component.

``--interaction``
   Interaction random effect ``<a>:<b>``, the rescaled Hadamard product of two
   kernels. Each operand is a loaded effect name or a GRM prefix.

.. rubric:: Input and Covariate Options

``--pheno-col`` ``0``
   0-based trait column index after ``FID``/``IID`` (first trait = 0).

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Discrete covariate TSV in format ``FID IID factor1 ...``.

.. rubric:: Model Configuration

``--mode`` ``A``
   Wald test mode: ``A`` (additive, single, df=1), ``D`` (dominance, single,
   df=1), or ``AD`` (joint additive + dominance, df=2).

``--geno-method, --gm`` ``OCH``
   Genotype coding method. Accepts codes ``SH``, ``CH``, ``OSH``, ``OCH``,
   ``S``, ``C``, ``OS``, ``OC``, ``NS``, ``NC``. See
   :ref:`genotype-processor-methods` for what each code means.

``--transform`` ``none``
   Phenotype transform: ``none``, ``dint`` (Direct INT), ``iint`` (Indirect INT).

``--int-offset`` ``0.375``
   Rank-INT offset parameter ``k``.

``--loco`` ``false``
   Use leave-one-chromosome-out GRMs. Cannot be combined with ``--interaction``.

.. rubric:: REML and Performance

``--max-iter`` ``100``
   Maximum model-fit (REML) iterations.

``--tol`` ``1e-06``
   Convergence tolerance.

``-c, --chunk-size`` ``10000``
   Number of SNPs per association-testing chunk.

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads to use.

Output Files
------------

After a successful run, GWAS summary statistics are written to:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File pattern
     - Contents
     - Reference
   * - ``<out>.gwas.tsv``
     - SNP-wise test statistics (effect size, SE, P-value, allele fields)
     - :ref:`gwas-output-format`
   * - ``<out>.log``
     - Run log with configuration, dataset summary, and timing.
     -

Warnings and Notes
------------------

.. warning::

   ``--loco`` requires chromosome-wise GRM inputs generated from
   ``gelex grm --loco``. Use the matching GRM prefix in ``--grm``. ``--loco``
   cannot be combined with ``--interaction``.

.. note::

   For ``--mode D`` or ``--mode AD``, provide GRM inputs consistent with the
   dominance model setup (typically additive + dominance GRMs).

Examples
--------

.. code-block:: bash
   :caption: Standard Additive GWAS

   gelex assoc \
      -b genotypes_qc \
      -p phenotypes.tsv \
      --grm my_grm \
      -o basic_gwas

.. code-block:: bash
   :caption: Add Quantitative and Categorical Covariates

   gelex assoc \
      -b genotypes_qc \
      -p phenotypes.tsv \
      --grm my_grm \
      --qcovar pcs.tsv \
      --dcovar sex.tsv \
      -o covar_gwas

.. code-block:: bash
   :caption: LOCO Analysis

   gelex assoc \
      -b genotypes_qc \
      -p phenotypes.tsv \
      --grm my_grm.add \
      --loco \
      -o loco_gwas

.. code-block:: bash
   :caption: Dominance Model with Two GRMs

   gelex assoc \
      -b genotypes_qc \
      -p phenotypes.tsv \
      --grm my_grm.add my_grm.dom \
      --mode D \
      --transform iint \
      -o dom_gwas

.. code-block:: bash
   :caption: Larger Thread Count

   gelex assoc \
      -b genotypes_qc \
      -p phenotypes.tsv \
      --grm my_grm \
      --threads 16 \
      -o fast_gwas

See Also
--------

- :ref:`grm-command` for preparing GRM inputs.
- :ref:`gwas-output-format` for GWAS output columns.
