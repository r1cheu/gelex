GWAS Tutorial
=============

.. admonition:: Quick Start
   :class: tip

   If you are in a hurry, here is the recommended "best practice" pipeline using `plink2` for preprocessing and **Gelex** for MLM-based association testing:

   .. code-block:: bash
      :caption: Quality Control (QC): Filter by MAF > 0.01 and Missingness < 0.05

      plink2 --bfile raw_data --maf 0.01 --geno 0.05 --make-bed --out genotypes_qc

   .. code-block:: bash
      :caption: LD Pruning (Calculate): Identify independent SNPs (window 50, step 10, r^2 0.1)

      plink2 --bfile genotypes_qc --indep-pairwise 50 10 0.1 --out pruned

   .. code-block:: bash
      :caption: LD Pruning (Extract): Create pruned dataset

      plink2 --bfile genotypes_qc --extract pruned.prune.in --make-bed --out genotypes_pruned

   .. code-block:: bash
      :caption: Population Structure: Compute Principal Components (PCs)

      plink2 --bfile genotypes_pruned --pca 20 --out population_pcs

   .. code-block:: bash
      :caption: Compute GRM: Use pruned SNPs for population structure

      gelex grm -b genotypes_pruned --add -o my_grm

   .. code-block:: bash
      :caption: Run Association Test: Test QC'ed SNPs controlling for structure

      gelex assoc \
        -b genotypes_qc \
        -p phenotypes.tsv \
        --grm my_grm \
        --qcovar population_pcs.eigenvec \
        --transform iint \
        -o final_gwas_results


Background
----------


Genome-wide association study (GWAS) tests the statistical association between genetic variants (SNPs) and a phenotype of interest across the entire genome. The goal is to identify SNPs that are significantly associated with the trait, which may point to causal genes or regulatory regions.

.. seealso::
   See :doc:`data_formats` for details on the files used in this process.


Mixed Linear Model (MLM)
~~~~~~~~~~~~~~~~~~~~~~~~

Gelex uses a **mixed linear model (MLM)** to account for population structure and cryptic relatedness:

.. math::

   \mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{x}b + \mathbf{g} + \boldsymbol{\varepsilon}

where:

- :math:`\mathbf{y}` is the phenotype vector.
- :math:`\mathbf{X}\boldsymbol{\beta}` represents fixed covariate effects.
- :math:`\mathbf{x}b` is the fixed effect of the candidate SNP.
- :math:`\mathbf{g}` is the random polygenic effect, :math:`\mathbf{g} \sim \mathcal{N}(\mathbf{0},\, \sigma_g^2 \mathbf{G})`.
- :math:`\boldsymbol{\varepsilon}` is the residual error.

**Efficiency:**
Similar to tools like GCTA, EMMAX, and GEMMA, Gelex estimates the variance components (:math:`\sigma_g^2` and :math:`\sigma_e^2`) once using a **Null Model** (excluding the candidate SNP). These estimates are then fixed when testing each SNP, drastically improving computational speed.

LOCO Strategy
~~~~~~~~~~~~~

The **Leave-One-Chromosome-Out (LOCO)** method excludes markers on the current chromosome :math:`k` from the GRM when testing SNPs on chromosome :math:`k`. This avoids proximal contamination and increases power.

**Implementation in Gelex:**

Gelex requires two sets of GRMs for LOCO analysis:
1. **Global GRM:** Calculated using all SNPs across the genome.
2. **Chromosome-specific GRMs:** Calculated using SNPs from each chromosome individually.

During analysis (``--loco``), Gelex loads the global GRM and dynamically subtracts the contribution of the specific chromosome GRM to compute the "LOCO GRM" (:math:`G_{-k}`).

.. math::
   G_{-k} = \frac{G_{whole} \cdot K_{whole} - G_k \cdot K_k}{K_{whole} - K_k}

Phenotype Transformation
~~~~~~~~~~~~~~~~~~~~~~~~
~

GWAS assumes approximately normally distributed residuals. For phenotypes that deviate from normality, Gelex provides **Inverse Normal Transformation (INT)** options:

- **Direct INT** (``dint``): Rank-based transformation applied directly to the raw phenotype values.
- **Indirect INT** (``iint``): Transformation applied to the residuals after regressing out covariate effects, preserving covariate-phenotype relationships.

Both methods use the Blom offset :math:`k = 3/8` by default: :math:`\Phi^{-1}\!\left(\frac{r_i - k}{n - 2k + 1}\right)`, where :math:`r_i` is the rank of individual :math:`i`.

Workflow Overview
~~~~~~~~~~~~~~~~~

A typical GWAS analysis in Gelex involves two steps:

1. **Compute GRM** — Build a genomic relationship matrix (with optional LOCO mode) using ``gelex grm``.
2. **Run association test** — Perform per-SNP association testing using ``gelex assoc`` with the pre-computed GRM.

Step 1: Compute GRM
--------------------

Before running association analysis, you need to compute a genomic relationship matrix (GRM) using the ``grm`` subcommand. For a full list of options, see :ref:`grm-command`. For more information on the GRM file structure, see :ref:`grm-format`.

Data Preparation
~~~~~~~~~~~~~~~~

.. tip::
   For computing the GRM and Principal Components (PCs), it is highly recommended to follow a **two-stage quality control process**:

   1. **Basic QC**: Filter out variants with low minor allele frequency (MAF) or high missingness.
   2. **LD Pruning**: Use a subset of SNPs that are in approximate Linkage Equilibrium (LD). This significantly reduces computational burden and avoids bias caused by regions of high LD.

We recommend using `plink2 <https://www.cog-genomics.org/plink/2.0/>`_ for these steps.

.. warning::
   **Sample ID Consistency**: Ensure that Sample IDs in your PLINK binary files match exactly with those in your phenotype and covariate files. Inconsistent IDs will result in samples being excluded from the analysis.

.. code-block:: bash
   :caption: Quality Control: Filter by MAF > 0.01 and Missingness < 0.05

   plink2 --bfile raw_data --maf 0.01 --geno 0.05 --make-bed --out genotypes_qc

.. code-block:: bash
   :caption: Identify SNPs to prune (window 50, step 10, r^2 0.1)

   plink2 --bfile genotypes_qc --indep-pairwise 50 10 0.1 --out ld_prune

.. code-block:: bash
   :caption: Create a pruned bfile for GRM and PCA

   plink2 --bfile genotypes_qc --extract ld_prune.prune.in --make-bed --out genotypes_pruned

.. code-block:: bash
   :caption: Compute Principal Components (PCs) for population stratification

   plink2 --bfile genotypes_pruned --pca 20 --out population_pcs

Basic Usage
~~~~~~~~~~~

.. code-block:: bash
   :caption: Compute GRM: Use pruned SNPs

   gelex grm -b genotypes_pruned --add -o my_grm

Step 2: Association Analysis
----------------------------

Run association testing using the ``assoc`` subcommand.

.. tip::
   **Which SNPs to use?**
   It is standard practice to test **all QC'ed SNPs** (``genotypes_qc``) while using a GRM built from **pruned SNPs** (``genotypes_pruned``) to control for population structure.

Examples
~~~~~~~~

**Standard MLM (Basic)**

This performs a standard mixed model association test.

.. code-block:: bash
   :caption: Standard MLM

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm \
     --qcovar population_pcs.eigenvec \
     -o basic_gwas

**LOCO Analysis**

LOCO analysis requires pre-computing both global and chromosome-specific GRMs.

.. code-block:: bash
   :caption: Compute Global GRM (Result: my_grm.add.bin)

   gelex grm -b genotypes_pruned --add -o my_grm

.. code-block:: bash
   :caption: Compute Chromosome GRMs (Result: my_grm.add.chr*.bin)

   gelex grm -b genotypes_pruned --add --loco -o my_grm

.. code-block:: bash
   :caption: Run Association with LOCO

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm.add \
     --qcovar population_pcs.eigenvec \
     --loco \
     -o loco_gwas

**With Covariates & Transformation**

Include discrete/quantitative covariates and apply Inverse Normal Transformation (INT) to residuals.

.. code-block:: bash
   :caption: Covariates & INT (using standard MLM)

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_grm \
     --qcovar population_pcs.eigenvec \
     --dcovar sex.tsv \
     --transform iint \
     -o refined_gwas

**Dominance Model**

Test for dominance effects using both additive and dominance GRMs.

.. code-block:: bash
   :caption: Compute additive & dominance GRMs using pruned SNPs

   gelex grm -b genotypes_pruned --add --dom -o my_background

.. code-block:: bash
   :caption: Run dominance association test using both GRMs

   gelex assoc \
     -b genotypes_qc \
     -p phenotypes.tsv \
     --grm my_background.add my_background.dom \
     --model d \
     --transform iint \
     -o dom_gwas

**High Performance (Multi-threading)**

For large datasets, you can increase the number of threads.

.. note::
   For optimal performance with dense linear algebra operations (MKL/Eigen), set ``--threads`` to the number of **physical cores**, not logical threads (hyperthreading). Using hyperthreads often degrades performance due to resource contention.

.. code-block:: bash
   :caption: Multi-threaded Analysis (e.g., for a 10-core CPU)

   gelex assoc -b genotypes_qc -p phenotypes.tsv --grm my_grm --threads 10 -o fast_gwas

Output Format
-------------

The association analysis produces a tab-separated file ``<out>.gwas.tsv``. See :ref:`gwas-output-format` for a detailed description of each column.

**Key Columns:**

BETA
  Effect size of the allele A1.

SE
  Standard error of the effect size.

P
  P-value from the Wald test.

**Example Output:**

.. code-block:: text

   CHR    SNP               BP      A1    A2    A1FREQ      BETA          SE            P
   1      chr01_1266_G_A    1266    A     G     0.335366    -0.0104849    0.0448972     8.153484e-01
   1      chr01_1325_C_T    1325    T     C     0.337159     0.0106326    0.0521935     8.385760e-01
   1      chr01_1335_G_T    1335    T     G     0.334469     0.0153185    0.0474137     7.466333e-01
   ...    ...               ...     ...   ...   ...         ...           ...           ...
