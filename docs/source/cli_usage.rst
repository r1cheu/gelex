CLI Usage
=========

Gelex provides a comprehensive command-line interface for genomic prediction and simulation. This documentation covers all available subcommands, options, and usage patterns.

Overview
--------

Gelex supports two main subcommands:

- ``fit`` - Fit genomic prediction models using Bayesian or GBLUP methods
- ``simulate`` - Simulate phenotypes based on genetic data

Basic Usage
~~~~~~~~~~~

.. code-block:: bash

   # Display help
   gelex --help

   # Fit command help
   gelex fit --help

   # Simulate command help
   gelex simulate --help

fit Command
-----------

The ``fit`` command performs genomic prediction using various Bayesian and frequentist methods.

Basic Syntax
~~~~~~~~~~~~

.. code-block:: bash

   gelex fit [OPTIONS]

Options
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Option
     - Default
     - Description / Format
   * - ``-p, --pheno``
     - (required)
     - Phenotype file (TSV: FID, IID, trait1, trait2, ...)
   * - ``-b, --bfile``
     - (required)
     - PLINK binary file prefix (requires .bed, .bim, .fam files)
   * - ``--pheno-col``
     - 2
     - Phenotype column index (0-based)
   * - ``--chunk-size``
     - 10000
     - SNPs per chunk (controls memory usage)
   * - ``--iid-only``
     - false
     - Use only IID for sample matching (ignore FID)
   * - ``-m, --method``
     - RR
     - Prediction method (see Method Reference below)
   * - ``--scale``
     - 0,0.001,0.01,0.1,1
     - Additive variance scales for BayesR
   * - ``--pi``
     - 0.95,0.05
     - Additive mixture proportions for BayesB/C/R
   * - ``--dscale``
     - 0,0.001,0.01,0.1,1
     - Dominance variance scales for BayesR
   * - ``--dpi``
     - 0.95,0.05
     - Dominance mixture proportions for BayesB/C/R
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
   * - ``--threads``
     - CPU cores / 2
     - Number of CPU threads to use
   * - ``--mmap``
     - false
     - Use memory-mapped I/O for genotype matrix
   * - ``--qcovar``
     - ""
     - Quantitative covariates (TSV: FID, IID, covar1, ...)
   * - ``--covar``
     - ""
     - Categorical covariates (TSV: FID, IID, factor1, ...)
   * - ``-o, --out``
     - "gelex"
     - Output file prefix

Method Reference
~~~~~~~~~~~~~~~~

Gelex supports the following methods for the ``--method`` option:

Bayesian Models (Additive Effects Only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Additive Bayesian Models
   :header-rows: 1
   :widths: 20 80

   * - Method
     - Description
   * - ``A``
     - BayesA: t-distributed marker effects
   * - ``B``
     - BayesB: Mixture model with spike-slab prior
   * - ``Bpi``
     - BayesB with mixture proportion estimation
   * - ``C``
     - BayesC: Similar to BayesB with different prior
   * - ``Cpi``
     - BayesC with mixture proportion estimation
   * - ``R``
     - BayesR: Multiple variance components
   * - ``RR``
     - BayesRR: Ridge regression (Bayesian LASSO)

Bayesian Models with Dominance Effects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Bayesian Models with Dominance
   :header-rows: 1
   :widths: 20 80

   * - Method
     - Description
   * - ``Ad``
     - BayesA with dominance effects
   * - ``Bd``
     - BayesB with dominance effects
   * - ``Bdpi``
     - BayesB with dominance and mixture estimation
   * - ``Cd``
     - BayesC with dominance effects
   * - ``Cdpi``
     - BayesC with dominance and mixture estimation
   * - ``Rd``
     - BayesR with dominance effects
   * - ``RRd``
     - BayesRR with dominance effects

Frequentist Models
^^^^^^^^^^^^^^^^^^

.. list-table:: Frequentist Models
   :header-rows: 1
   :widths: 20 80

   * - Method
     - Description
   * - ``GBLUP``
     - Genomic Best Linear Unbiased Prediction

For detailed information about models, see the :doc:`bayesian_models` and :doc:`frequentist_models` documentation.

simulate Command
----------------

The ``simulate`` command generates simulated phenotypes based on genetic data and specified parameters.

Basic Syntax
~~~~~~~~~~~~

.. code-block:: bash

   gelex simulate [OPTIONS]

Options
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Option
     - Description
     - Format
   * - ``-b, --bfile <BFILE>``
     - PLINK binary file prefix
     - Requires .bed, .bim, .fam files with same prefix
   * - ``-c, --causal <CAUSAL>``
     - Causal variants file
     - TSV: SNP ID per line, optional effect size
   * - ``--h2 <double>``
     - 0.5
     - Narrow-sense heritability (range: 0-1)
   * - ``--seed <int>``
     - -1
     - Random seed for reproducibility (-1 for time-based)
   * - ``-o, --out <OUT>``
     - "sim.phen"
     - Output file prefix for simulated phenotypes

Usage Examples
--------------

Basic Bayesian Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   gelex fit \
     --bfile genotypes \
     --pheno phenotypes.tsv \
     --method RR \
     --iters 10000 \
     --burnin 2000 \
     --threads 8 \
     --out my_analysis

Bayesian Model with Dominance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   gelex fit \
     --bfile genotypes \
     --pheno phenotypes.tsv \
     --method RRd \
     --iters 8000 \
     --burnin 3000 \
     --thin 5 \
     --chains 3 \
     --out dominance_analysis

With Covariates
~~~~~~~~~~~~~~~

.. code-block:: bash

   gelex fit \
     --bfile genotypes \
     --pheno phenotypes.tsv \
     --qcovar quantitative_covariates.tsv \
     --covar categorical_covariates.tsv \
     --method Bpi \
     --iters 12000 \
     --out analysis_with_covariates

Memory-Efficient Large Dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   gelex fit \
     --bfile large_genotypes \
     --pheno phenotypes.tsv \
     --method RR \
     --chunk-size 5000 \
     --mmap \
     --threads 4 \
     --out large_analysis

Phenotype Simulation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   gelex simulate \
     --bfile genotypes \
     --causal causal_variants.tsv \
     --h2 0.3 \
     --seed 12345 \
     --out simulated_phenotypes
