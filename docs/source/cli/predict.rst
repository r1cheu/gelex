.. _predict-command:

predict
=======

Generate genomic predictions using fitted SNP effects.

Basic Syntax
------------

.. code-block:: bash
   :caption: Basic Usage

   gelex predict --bfile <bfile> --snp-eff <file> --out <output> [OPTIONS]

Options
-------

.. rubric:: Data Files

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``-b, --bfile``
     - (required)
     - PLINK binary file prefix (.bed/.bim/.fam)
   * - ``-e, --snp-eff``
     - (required)
     - SNP effects file (.snp.eff)
   * - ``--covar-eff``
     - ""
     - Covariate effects file (.param)
   * - ``--qcovar``
     - ""
     - Quantitative covariates file
   * - ``--dcovar``
     - ""
     - Discrete covariates file
   * - ``-o, --out``
     - (required)
     - Output file path for predictions

.. rubric:: Processing Options

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Option
     - Default
     - Description
   * - ``--iid-only``
     - false
     - Use only IID for sample matching (ignore FID)
   * - ``-c, --chunk-size``
     - 10000
     - SNPs per chunk (controls memory usage)

Examples
--------

.. code-block:: bash
   :caption: Basic Genomic Prediction

   gelex predict \
     -b target_data \
     -e trained_model.snp.eff \
     -o predictions.tsv

.. code-block:: bash
   :caption: Prediction with Covariate Effects

   gelex predict \
     -b target_data \
     -e trained_model.snp.eff \
     --covar-eff trained_model.param \
     --qcovar target_age.tsv \
     --dcovar target_sex.tsv \
     -o predictions.tsv

.. code-block:: bash
   :caption: IID-Only Matching

   gelex predict \
     -b target_data \
     -e trained_model.snp.eff \
     --iid-only \
     -o predictions.tsv

.. code-block:: bash
   :caption: Low Memory Mode

   gelex predict \
     -b target_data \
     -e trained_model.snp.eff \
     --chunk-size 1000 \
     -o predictions.tsv
