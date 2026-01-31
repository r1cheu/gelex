.. _fit-command:

fit
===

The ``fit`` command performs genomic prediction using various Bayesian methods (BayesAlphabet).

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex fit --pheno <pheno_file> --bfile <genotype_prefix> --method <method> [OPTIONS]

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
   * - ``-b, --bfile``
     - (required)
     - PLINK binary file prefix (.bed/.bim/.fam)
   * - ``--qcovar``
     - ""
     - Quantitative covariates (TSV: FID, IID, covar1, ...)
   * - ``--dcovar``
     - ""
     - Discrete (categorical) covariates (TSV: FID, IID, factor1, ...)
   * - ``-o, --out``
     - "gelex"
     - Output file prefix

.. rubric:: Processing Options

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--pheno-col``
     - 2
     - Phenotype column index (0-based)
   * - ``-c, --chunk-size``
     - 10000
     - SNPs per chunk (controls memory usage)
   * - ``--iid-only``
     - false
     - Use only IID for sample matching (ignore FID)

.. rubric:: Model Configuration

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-m, --method``
     - RR
     - Method: A/B/C/R/RR (+d for dominance, +pi to estimate mixture)
   * - ``--scale``
     - 0,0.001,...
     - Additive variance scales for BayesR (5 values)
   * - ``--pi``
     - 0.95,0.05
     - Additive mixture proportions for BayesB/C/R
   * - ``--dscale``
     - 0,0.001,...
     - Dominance variance scales for BayesR (5 values)
   * - ``--dpi``
     - 0.95,0.05
     - Dominance mixture proportions for BayesB/C/R

.. rubric:: MCMC Configuration

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--iters``
     - 5000
     - Total MCMC iterations
   * - ``--burnin``
     - 4000
     - Burn-in iterations to discard
   * - ``--thin``
     - 1
     - Thinning interval for samples
   * - ``--chains``
     - 1
     - Number of MCMC chains

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
   * - ``--mmap``
     - false
     - Use memory-mapped I/O (much lower RAM, may be slower)

Examples
--------

.. code-block:: bash
   :caption: BayesRR (Ridge Regression) - Default

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m RR \
     -o model_rr

.. code-block:: bash
   :caption: BayesR (Mixture Model) - Recommended for accuracy

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m R \
     --chains 4 \
     -o model_bayesr

.. code-block:: bash
   :caption: BayesB (Variable Selection)

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m B \
     --pi 0.99 0.01 \
     -o model_bayesb

.. code-block:: bash
   :caption: With Covariates (Fixed Effects)

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m R \
     --dcovar sex.tsv \
     --qcovar age.tsv \
     -o model_covar

.. code-block:: bash
   :caption: High Precision Training

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m R \
     --iters 50000 \
     --burnin 10000 \
     --thin 5 \
     -o model_high_prec

.. code-block:: bash
   :caption: Dominance Model (BayesRd)

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m Rd \
     --dscale 0.0001 0.001 0.01 0.1 1.0 \
     --dpi 0.95 0.05 \
     -o model_dom

.. code-block:: bash
   :caption: Parameter Tuning (BayesCpi with Custom Priors)

   gelex fit \
     -b train_data \
     -p phenotypes.tsv \
     -m Cpi \
     --pi 0.9 0.1 \
     --scale 0.01 0.1 \
     -o model_cpi
