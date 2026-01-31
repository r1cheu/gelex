Genomic Selection Tutorial
==========================

.. admonition:: Quick Start
   :class: tip

   For a standard Genomic Selection pipeline using **BayesR**:

   .. code-block:: bash
      :caption: Model Fitting (Train)

      gelex fit \
        --bfile train_data \
        --pheno train_pheno.tsv \
        --method R \
        --iters 20000 --burnin 5000 \
        --out trained_model

   .. code-block:: bash
      :caption: Genomic Prediction (Predict)

      gelex predict \
        --bfile target_data \
        --snp-eff trained_model.snp.eff \
        --out predictions.tsv

Background
----------

Genomic Selection (GS) uses genome-wide markers to predict complex traits. Unlike GWAS, which focuses on identifying specific significant variants, GS aims to capture the total genetic value of an individual by simultaneously estimating the effects of all markers, even those with small effects.

Gelex implements the **BayesAlphabet** family of models (BayesA, BayesB, BayesC, BayesR, BayesRR), which use different prior distributions to model the genetic architecture of traits ranging from simple (oligogenic) to complex (polygenic).

Workflow Overview
-----------------

A typical GS analysis involves two main steps:

1.  **Model Fitting (`fit`)**: Train a Bayesian model on a reference population with both genotypes and phenotypes to estimate marker effects.
2.  **Prediction (`predict`)**: Apply the estimated marker effects to the genotypes of a target population (candidates) to predict their Genomic Estimated Breeding Values (GEBVs) or phenotypic values.

Step 1: Model Fitting
---------------------

The first step is to fit a model to your training data. This estimates the effect size of each SNP and any fixed covariates.

.. seealso::
   See :ref:`fit-command` for a full list of options.

Choose a Method
~~~~~~~~~~~~~~~

*   **BayesRR / Ridge Regression**: Assumes all SNPs have non-zero effects drawn from a normal distribution. Good for highly polygenic traits.
*   **BayesA**: Assumes all SNPs have non-zero effects but allows for different variances (heavy-tailed).
*   **BayesB / BayesC**: Variable selection models that assume a proportion of SNPs have zero effect (`pi` parameter).
*   **BayesR**: A flexible mixture model that assumes SNP effects come from a mixture of normal distributions with different variances (e.g., zero, small, medium, large). Often provides the highest accuracy across diverse genetic architectures.

Basic Usage
~~~~~~~~~~~

To fit a **BayesR** model:

.. code-block:: bash
   :caption: Fit BayesR Model

   gelex fit \
     --bfile train_genotypes \
     --pheno phenotypes.tsv \
     --method R \
     --chains 4 \
     --out model_output

Handling Covariates
~~~~~~~~~~~~~~~~~~~

You can include fixed effects such as sex (discrete) or age (quantitative):

.. code-block:: bash
   :caption: Fit Model with Covariates

   gelex fit \
     --bfile train_genotypes \
     --pheno phenotypes.tsv \
     --dcovar sex.tsv \
     --qcovar age.tsv \
     --method R \
     --out model_with_covars

Outputs
~~~~~~~

The `fit` command generates several files:

*   ``<out>.snp.eff``: Estimated SNP effects (used for prediction).
*   ``<out>.param``: Estimated hyper-parameters and covariate effects.
*   ``<out>.log``: Log of the MCMC process.

Step 2: Genomic Prediction
--------------------------

Once you have the estimated SNP effects (and optionally covariate effects), you can predict phenotypes for a target population.

.. seealso::
   See :ref:`predict-command` for detailed options.

Basic Usage
~~~~~~~~~~~

.. code-block:: bash
   :caption: Basic Prediction

   gelex predict \
     --bfile target_genotypes \
     --snp-eff model_output.snp.eff \
     --out predicted_values.tsv

Using Covariate Effects
~~~~~~~~~~~~~~~~~~~~~~~

If your training model included covariates, you can also use those estimated effects for prediction, provided the target population has the corresponding covariate data:

.. code-block:: bash
   :caption: Prediction with Covariates

   gelex predict \
     --bfile target_genotypes \
     --snp-eff model_with_covars.snp.eff \
     --covar-eff model_with_covars.param \
     --qcovar target_age.tsv \
     --dcovar target_sex.tsv \
     --out predicted_values_with_covars.tsv

Output Format
-------------

The output file is a TSV containing the predicted values:

.. code-block:: text
   :caption: Example Output

   FID    IID    PRS          Total_Value
   fam1   ind1   1.234e-01    1.543e-01
   fam1   ind2   -0.456e-01   -0.123e-01
   ...

*   **PRS**: Polygenic Risk Score (sum of SNP effects).
*   **Total_Value**: PRS + Covariate Effects (if provided).
