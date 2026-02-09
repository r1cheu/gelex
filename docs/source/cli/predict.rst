.. _predict-command:

predict
=======

Generate genomic predictions from fitted SNP effects.

Use this command after :ref:`fit-command` to score target samples.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex predict -b target_data -e model.snp.eff -o target.pred.tsv

.. code-block:: bash
   :caption: Full Syntax Template

   gelex predict --bfile <genotype_prefix> --snp-eff <snp_effect_file> --out <prediction_file> [OPTIONS]

Required inputs are genotype prefix (``--bfile``), SNP effects
(``--snp-eff``), and output path (``--out``).

Options
-------

.. rubric:: Quick Start Options

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``) for target samples.

``-e, --snp-eff`` ``required``
   SNP effects file from ``gelex fit`` (usually ``<out>.snp.eff``).

``-o, --out`` ``required``
   Output file path for prediction results.

.. rubric:: Input and Covariate Options

``--covar-eff``
   Optional covariate effect file (usually ``<out>.param`` from ``fit``).

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Categorical covariate TSV in format ``FID IID factor1 ...``.

.. rubric:: Processing Options

``--iid-only`` ``false``
   Match samples by IID only and ignore FID.

``-c, --chunk-size`` ``10000``
   Number of SNPs per processing chunk. Lower values reduce peak memory.

Output Files
------------

``predict`` writes one prediction table at the exact path passed to ``--out``.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File
     - Contents
     - Reference
   * - ``<out>``
     - Per-sample predictions (FID/IID and predicted values)
     - :ref:`predict-output-format`

Warnings and Notes
------------------

.. warning::

   Use ``--iid-only`` only when IID uniquely identifies individuals in every
   file. Otherwise, samples may be mismatched.

.. note::

   If you use ``--covar-eff``, keep covariate files consistent with the fit
   stage (same variables, coding, and compatible sample IDs).

Examples
--------

.. code-block:: bash
   :caption: Basic Genomic Prediction

   gelex predict \
      -b target_data \
      -e trained_model.snp.eff \
      -o predictions.pred.tsv

.. code-block:: bash
   :caption: Add Covariate Effects

   gelex predict \
      -b target_data \
      -e trained_model.snp.eff \
      --covar-eff trained_model.param \
      --qcovar target_age.tsv \
      --dcovar target_sex.tsv \
      -o predictions_with_covar.pred.tsv

.. code-block:: bash
   :caption: IID-Only Matching

   gelex predict \
      -b target_data \
      -e trained_model.snp.eff \
      --iid-only \
      -o predictions_iid_only.pred.tsv

.. code-block:: bash
   :caption: Low-Memory Chunking

   gelex predict \
      -b target_data \
      -e trained_model.snp.eff \
      --chunk-size 2000 \
      -o predictions_low_mem.pred.tsv

See Also
--------

- :ref:`fit-command` for training SNP and covariate effects.
- :ref:`predict-output-format` for prediction output columns.
