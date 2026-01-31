.. _assoc-command:

assoc
=====

Perform genome-wide association study (GWAS) using mixed linear model.

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex assoc --pheno <pheno> --bfile <bfile> --grm <grm_prefix> [OPTIONS]

Options
-------

.. rubric:: Data Files

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-p, --pheno``
     - (required)
     - Phenotype file (TSV format: FID, IID, trait1, ...)
   * - ``--pheno-col``
     - 2
     - Phenotype column index (0-based)
   * - ``-b, --bfile``
     - (required)
     - PLINK binary file prefix (.bed/.bim/.fam)
   * - ``--grm``
     - (required)
     - GRM file prefix(es). Can specify multiple GRMs.
   * - ``--qcovar``
     - ""
     - Quantitative covariates (TSV: FID, IID, covar1, ...)
   * - ``--dcovar``
     - ""
     - Discrete (categorical) covariates (TSV: FID, IID, factor1, ...)
   * - ``-o, --out``
     - "gelex"
     - Output file prefix

.. rubric:: REML Options

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--max-iter``
     - 100
     - Max iteration in REML process
   * - ``--tol``
     - 1e-06
     - tolerance for convergence in REML process

.. rubric:: Processing Options

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-c, --chunk-size``
     - 10000
     - SNPs per chunk for association testing
   * - ``--iid-only``
     - false
     - Use only IID for sample matching (ignore FID)
   * - ``--loco``
     - false
     - Enable Leave-One-Chromosome-Out (LOCO) mode

.. rubric:: Model Configuration

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--model``
     - a
     - Association model: ``a`` (additive), ``d`` (dominance)
   * - ``--transform``
     - none
     - Phenotype transformation: ``none``, ``dint``, ``iint``
   * - ``--int-offset``
     - 0.375
     - INT offset parameter k (default: 3/8 Blom offset)

.. rubric:: Performance

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-t, --threads``
     - 12
     - Number of CPU threads to use

Examples
--------

.. code-block:: bash
   :caption: Standard MLM (Basic)

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm \
     --qcovar population_pcs.eigenvec \
     -o basic_gwas

.. code-block:: bash
   :caption: LOCO Analysis

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm.add \
     --qcovar population_pcs.eigenvec \
     --loco \
     -o loco_gwas

.. code-block:: bash
   :caption: With Covariates & Transformation (INT)

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm \
     --qcovar population_pcs.eigenvec \
     --dcovar sex.tsv \
     --transform iint \
     -o refined_gwas

.. code-block:: bash
   :caption: Dominance Model

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm.add my_grm.dom \
     --model d \
     --transform iint \
     -o dom_gwas

.. code-block:: bash
   :caption: High Performance (Multi-threading)

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm \
     --threads 16 \
     -o fast_gwas
