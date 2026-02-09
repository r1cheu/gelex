.. _fit-command:

fit
===

Train SNP effect models for genomic prediction using BayesAlphabet methods.

Use this command when you want to learn marker effects from training data,
then reuse those effects in :ref:`predict-command`.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex fit -b train_data -p phenotypes.tsv -m RR -o model_rr

.. code-block:: bash
   :caption: Full Syntax Template

   gelex fit --pheno <pheno_file> --bfile <genotype_prefix> --method <method> [OPTIONS]

Required inputs are phenotype file (``--pheno``), genotype prefix (``--bfile``),
and model method (``--method``).

Method Selection
----------------

Choose a method based on your goal before tuning other parameters.

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Method
     - Use when
     - Trade-off
   * - ``RR``
     - All SNPs are assumed to have non-zero effects; use as a baseline.
     - Stable and simple, but weak variable selection.
   * - ``R``
     - You expect a mixture of effect sizes and want flexible shrinkage.
     - Better accuracy in many traits, with moderate runtime.
   * - ``B`` / ``C``
     - You expect many near-zero SNP effects and want explicit variable selection.
     - Stronger sparsity, but more sensitive to prior settings.
   * - ``A``
     - You want all SNPs included with SNP-specific shrinkage.
     - More MCMC sampling cost than ``RR``.
   * - ``Rd``
     - You want to model dominance effects alongside additive effects.
     - More parameters and longer runtime.
   * - ``Bpi`` / ``Cpi`` / ``Rpi``
     - You want the model to estimate mixture proportions from the data.
     - More adaptive, but may require longer chains for stable estimates.

If you are unsure, start with ``RR`` to establish a baseline, then try
``R`` as a stronger default for production runs.

Options
-------

.. rubric:: Quick Start Options

``-p, --pheno`` ``required``
   Phenotype TSV file (``FID IID trait1 ...``).

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-m, --method`` ``RR``
   Modeling method. Start with ``RR`` (baseline) or ``R``
   (accuracy-oriented).

``-o, --out`` ``gelex``
   Output prefix for generated files.

.. rubric:: Input Options

``-p, --pheno`` ``required``
   Phenotype TSV file in format ``FID IID trait1 ...``.

``--pheno-col`` ``2``
   0-based trait column index in the phenotype file.

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Categorical covariate TSV in format ``FID IID factor1 ...``.

``--iid-only`` ``false``
   Match samples by IID only and ignore FID.

.. rubric:: Model Options

``-m, --method`` ``RR``
   BayesAlphabet method. Supported: ``A/B/C/R/RR``; add ``d`` for dominance
   (for example ``Rd``); add ``pi`` to estimate mixture proportions (for
   example ``Cpi``).

``--scale`` ``0 0.001 0.01 0.1 1``
   Additive variance scales, typically used in BayesR-style models.

``--pi`` ``0.99 0.01``
   Additive mixture proportions. For BayesR, default is
   ``0.99 0.005 0.003 0.001 0.001``.

``--dscale`` ``0 0.001 0.01 0.1 1``
   Dominance variance scales for dominance-enabled models.

``--dpi`` ``0.99 0.01``
   Dominance mixture proportions. For BayesR dominance models, default is
   ``0.99 0.005 0.003 0.001 0.001``.

.. rubric:: MCMC Options

``--iters`` ``5000``
   Total MCMC iterations.

``--burnin`` ``4000``
   Initial iterations discarded before sampling.

``--thin`` ``1``
   Keep one sample every ``thin`` iterations.

``--chains`` ``1``
   Number of independent MCMC chains.

.. rubric:: Performance and Output

``-c, --chunk-size`` ``10000``
   Number of SNPs per processing chunk. Lower values reduce peak memory.

``-t, --threads`` ``12``
   Number of CPU threads to use.

``--mmap`` ``false``
   Enable memory-mapped I/O. Usually lowers RAM pressure and may reduce speed.

``-o, --out`` ``gelex``
   Output prefix for all generated files.

Output Files
------------

After a successful run, check files with your output prefix first.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File pattern
     - Contents
     - Typical next step
   * - ``<out>.snp.eff``
     - Estimated SNP effects
     - Use with ``gelex predict --snp-eff``
   * - ``<out>.param``
     - Estimated fixed/covariate effects and model parameters
     - Optional input for ``gelex predict --covar-eff``
   * - ``<out>*``
     - Run logs and model-specific artifacts
     - Review convergence and configuration used

Warnings and Notes
------------------

.. note::

   For many datasets, a practical starting point is ``--burnin`` around
   20%-50% of ``--iters``. Increase ``--iters`` when posterior summaries are
   unstable across runs.

.. warning::

   Use ``--iid-only`` only when IID uniquely identifies individuals in all
   files. If FID+IID pairs are required for your dataset, enabling this flag
   can silently mismatch samples.

.. note::

   If memory is limited, reduce ``--chunk-size`` first, then enable
   ``--mmap``. This usually lowers RAM usage with a possible runtime penalty.

Examples
--------

.. code-block:: bash
   :caption: Quick Start Baseline (RR)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m RR \
      -o model_rr

Expected outputs: ``model_rr.snp.eff``, ``model_rr.param``.

.. code-block:: bash
   :caption: Accuracy-Oriented Training (R)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --chains 4 \
      -o model_bayesr

Expected outputs: ``model_bayesr.snp.eff``, ``model_bayesr.param``.

.. code-block:: bash
   :caption: Sparse Effects with Variable Selection (B)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m B \
      --pi 0.99 0.01 \
      -o model_bayesb

.. code-block:: bash
   :caption: Add Fixed Effects (qcovar + dcovar)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --dcovar sex.tsv \
      --qcovar age.tsv \
      -o model_covar

.. code-block:: bash
   :caption: Longer MCMC for Stable Posterior Estimates

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --iters 50000 \
      --burnin 10000 \
      --thin 5 \
      -o model_high_prec

.. code-block:: bash
   :caption: Additive + Dominance Model (Rd)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m Rd \
      --dscale 0.0001 0.001 0.01 0.1 1.0 \
      --dpi 0.95 0.05 \
      -o model_dom

.. code-block:: bash
   :caption: Estimate Mixture Proportions (Cpi)

   gelex fit \
      -b train_data \
      -p phenotypes.tsv \
      -m Cpi \
      --pi 0.9 0.1 \
      --scale 0.01 0.1 \
      -o model_cpi

See Also
--------

- :ref:`predict-command` for applying trained effects to target samples.
- :ref:`assoc-command` for SNP-wise association analysis.
- :ref:`grm-command` for genomic relationship matrix construction.
