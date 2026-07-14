.. _predict-command:

predict
=======

Generate genomic predictions from fitted SNP effects.

Use this command after :ref:`mcmc-command` to score target samples.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex predict -b target_data -g trained_model -o target_predictions

.. code-block:: bash
   :caption: Full Syntax Template

   gelex predict --bfile <genotype_prefix> --gfile <model_prefix> --out <output_prefix> [OPTIONS]

Required inputs are the genotype prefix (``--bfile``), the fitted-model
prefix (``--gfile``), and the output prefix (``--out``).

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``) for target samples.

``-g, --gfile`` ``required``
   Fitted-model prefix from ``gelex mcmc`` or ``gelex vi``. Reads
   ``<prefix>.snpeff``, ``<prefix>.snpstats``, and ``<prefix>.param``
   (covariate coefficients are taken from ``.param``).

``-o, --out`` ``required``
   Output prefix. ``predict`` writes ``<out>.pred.tsv``.

.. rubric:: Covariate Options

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Categorical covariate TSV in format ``FID IID factor1 ...``.

.. rubric:: Processing Options

``-c, --chunk-size`` ``10000``
   Number of SNPs per processing chunk. Lower values reduce peak memory.

Output Files
------------

``predict`` writes one prediction table, using ``--out`` as a prefix.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File
     - Contents
     - Reference
   * - ``<out>.pred.tsv``
     - Per-sample predictions (FID, IID, prediction, plus covariate and
       per-mode GEBV columns)
     - :ref:`predict-output-format`

Warnings and Notes
------------------

.. note::

   Covariate coefficients are read from the model's ``.param`` file. When you
   pass ``--qcovar``/``--dcovar``, keep the covariate files consistent with the
   fit stage (same variables, coding, and compatible sample IDs).

Examples
--------

.. code-block:: bash
   :caption: Basic Genomic Prediction

   gelex predict \
      -b target_data \
      -g trained_model \
      -o predictions

.. code-block:: bash
   :caption: With Covariates

   gelex predict \
      -b target_data \
      -g trained_model \
      --qcovar target_age.tsv \
      --dcovar target_sex.tsv \
      -o predictions_with_covar

.. code-block:: bash
   :caption: Low-Memory Chunking

   gelex predict \
      -b target_data \
      -g trained_model \
      --chunk-size 2000 \
      -o predictions_low_mem

See Also
--------

- :ref:`mcmc-command` for training SNP and covariate effects.
- :ref:`predict-output-format` for prediction output columns.
