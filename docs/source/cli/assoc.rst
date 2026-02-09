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

Required inputs are phenotype file (``--pheno``), genotype prefix (``--bfile``),
and at least one GRM prefix (``--grm``).

Method Selection
----------------

Pick model and preprocessing strategy before tuning runtime options.

.. list-table::
   :header-rows: 1
   :widths: 22 43 35

   * - Option
     - Use when
     - Trade-off
   * - ``--model a``
      - You are running a standard additive-effect GWAS.
      - Fast and robust default for most analyses.
   * - ``--model d``
      - You want to test dominance effects in addition to additive effects.
      - Requires compatible GRMs and often larger sample sizes.
   * - ``--transform none``
      - The phenotype is already approximately normal.
      - Keeps interpretation on the original trait scale.
   * - ``--transform dint`` / ``iint``
      - The phenotype distribution is skewed or heavy-tailed.
      - Often improves calibration, but effect sizes are on transformed scale.
   * - ``--geno-method center``
      - You want the default genotype preprocessing pipeline.
      - Good default for stability and comparability across runs.

.. warning::

   If you use ``--transform`` with ``--model d``, dominance signals may be
   attenuated, which can reduce power to detect dominance effects.


Options
-------

.. rubric:: Quick Start Options

``-p, --pheno`` ``required``
   Phenotype TSV file in format ``FID IID trait1 ...``.

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``--grm`` ``required``
   One or more GRM prefixes.

``-o, --out`` ``gelex``
   Output prefix for GWAS results.

.. rubric:: Input and Covariate Options

``--pheno-col`` ``2``
   0-based trait column index in phenotype file.

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Categorical covariate TSV in format ``FID IID factor1 ...``.

``--iid-only`` ``false``
   Match samples by IID only and ignore FID.

.. rubric:: Model Configuration

``--model`` ``a``
   Association model: ``a`` (additive) or ``d`` (dominance).

``--geno-method`` ``center``
   Genotype method (center-family): ``center``, ``orth-center``,
   ``center-sample``, ``orth-center-sample``.

``--transform`` ``none``
   Phenotype transform: ``none``, ``dint`` (Direct INT), ``iint`` (Indirect INT).

``--int-offset`` ``0.375``
   INT offset parameter ``k``.

.. rubric:: REML and Performance

``--max-iter`` ``100``
   Maximum REML iterations.

``--tol`` ``1e-06``
   REML convergence tolerance.

``-c, --chunk-size`` ``10000``
   Number of SNPs per association-testing chunk.

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads to use.

``--loco`` ``false``
   Enable leave-one-chromosome-out analysis.

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

Warnings and Notes
------------------

.. warning::

   ``--loco`` requires chromosome-wise GRM inputs generated from
   ``gelex grm --loco``. Use the matching GRM prefix in ``--grm``.

.. note::

   For ``--model d``, provide GRM inputs consistent with the dominance model
   setup (typically additive + dominance GRMs).

.. warning::

   Use ``--iid-only`` only when IID uniquely identifies all individuals.

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
      --model d \
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
