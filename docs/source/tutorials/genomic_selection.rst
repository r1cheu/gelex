Genomic Selection
=================

.. admonition:: Quick Start
   :class: tip

   For a standard Genomic Selection pipeline using **BayesR**:

   .. code-block:: bash
      :caption: Model Fitting (Train)

      gelex mcmc \
        --bfile train_data \
        --pheno train_pheno.tsv \
        --method R \
        --iters 20000 --burn-in 5000 \
        --out trained_model

   .. code-block:: bash
      :caption: Genomic Prediction (Predict)

      gelex predict \
        --bfile target_data \
        --gfile trained_model \
        --out predictions

.. seealso::
   New to genomic selection or the BayesAlphabet models? See
   :doc:`/concepts/bayesian_models` for the background and how the methods differ.

Workflow Overview
-----------------

A typical GS analysis involves two main steps:

1.  **Model Fitting (``mcmc``)**: Train a Bayesian model on a reference population with both genotypes and phenotypes to estimate marker effects.
2.  **Prediction (``predict``)**: Apply the estimated marker effects to the genotypes of a target population (candidates) to predict their Genomic Estimated Breeding Values (GEBVs) or phenotypic values.

Step 1: Model Fitting
---------------------

The first step is to fit a model to your training data. This estimates the effect size of each SNP and any fixed covariates.

.. seealso::
   See :ref:`mcmc-command` for a full list of options.

Choose a Method
~~~~~~~~~~~~~~~

Select the model with ``--method`` (``RR``, ``A``, ``B``, ``C``, ``R``, ``CD``). Each
encodes a different assumption about the trait's genetic architecture; see
:doc:`/concepts/bayesian_models` for a comparison. When unsure, use
**BayesR** (``--method R``).

Basic Usage
~~~~~~~~~~~

To fit a **BayesR** model:

.. code-block:: bash
   :caption: Fit BayesR Model

   gelex mcmc \
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

   gelex mcmc \
     --bfile train_genotypes \
     --pheno phenotypes.tsv \
     --dcovar sex.tsv \
     --qcovar age.tsv \
     --method R \
     --out model_with_covars

Outputs
~~~~~~~

The ``mcmc`` command generates several files sharing the ``--out`` prefix:

*   ``<out>.snpeff``: Estimated SNP effects (used for prediction).
*   ``<out>.snpstats``: Per-SNP genotype statistics (used for prediction).
*   ``<out>.param``: Estimated hyper-parameters and covariate effects.
*   ``<out>.summary``: Posterior summary of model parameters.
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
     --gfile model_output \
     --out predicted_values

.. note::
   ``--gfile`` takes the fitted-model **prefix** (the ``mcmc`` ``--out`` value),
   not a single file. ``predict`` reads ``<prefix>.snpeff``, ``<prefix>.snpstats``
   and ``<prefix>.param`` from it. ``--out`` is likewise a prefix: the result is
   written to ``<out>.pred.tsv``.

Using Covariate Effects
~~~~~~~~~~~~~~~~~~~~~~~

If your training model included covariates, their estimated effects are stored in
``<prefix>.param`` and are applied automatically. You only need to supply the
target population's covariate data so it can be matched to those effects:

.. code-block:: bash
   :caption: Prediction with Covariates

   gelex predict \
     --bfile target_genotypes \
     --gfile model_with_covars \
     --qcovar target_age.tsv \
     --dcovar target_sex.tsv \
     --out predicted_values_with_covars

Output Format
-------------

The prediction is written to ``<out>.pred.tsv``, a tab-separated file with one row
per target individual:

.. code-block:: text
   :caption: Example Output (additive model)

   FID    IID    prediction    A
   fam1   ind1   0.123456      0.123456
   fam1   ind2   -0.045600     -0.045600
   ...

*   **FID** / **IID**: Family and individual identifiers.
*   **prediction**: Total predicted value (genetic value plus any covariate effects).
*   **A** / **D**: Genomic estimated value from the additive and, when the model
    includes dominance, dominance components.

When covariates are supplied, one column per covariate is inserted between
``prediction`` and the genetic-value columns, reporting that covariate's
contribution to the prediction.
