.. _mcmc-command:

mcmc
====

Train SNP effect models for genomic prediction using BayesAlphabet methods
with MCMC (Gibbs sampling).

Use this command when you want to learn marker effects from training data,
then reuse those effects in :ref:`predict-command`.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex mcmc -b train_data -p phenotypes.tsv -m RR -o model_rr

.. code-block:: bash
   :caption: Full Syntax Template

   gelex mcmc --pheno <pheno_file> --bfile <genotype_prefix> --method <method> [OPTIONS]

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
   * - ``CD``
     - You want coupled additive and dominance marker allocation.
     - More parameters and longer runtime.
   * - ``--sample-pi`` / ``--sample-dpi`` / ``--sample-jpi``
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

.. rubric:: Model Options

``-m, --method`` ``RR``
   BayesAlphabet method. Supported: ``RR``, ``A``, ``B``, ``C``, ``R``, ``CD``.

``--mode`` ``A``
   Genetic effect mode: ``A`` (additive), ``D`` (dominance), or ``AD``.
   ``CD`` requires ``AD``.

``--geno-method`` ``OrthStandardizeHWE``
   Genotype processing method. Available methods:
   ``StandardizeHWE`` (``SH``), ``CenterHWE`` (``CH``),
   ``OrthStandardizeHWE`` (``OSH``), ``OrthCenterHWE`` (``OCH``),
   ``Standardize`` (``S``), ``Center`` (``C``),
   ``OrthStandardize`` (``OS``), ``OrthCenter`` (``OC``).
   Abbreviations accepted.
   See :ref:`genotype-processor-methods`.

``--h2`` ``0.5``
   Additive heritability.

``--d2`` ``0.2``
   Dominance heritability.

``--random-pve``
   Variance fraction for non-SNP random effects.

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

``--jpi`` ``0.99 0.0033 0.0033 0.0033``
   Joint allocation proportions for ``CD``:
   both-off, additive-only, dominance-only, both-on.

``--sample-pi`` / ``--sample-dpi`` / ``--sample-jpi``
   Sample additive, dominance, or joint mixture proportions.

``--dom-pos-prob`` ``0.5``
   Initial probability that an active dominance effect is positive.

.. rubric:: MCMC Options

``--iters`` ``5000``
   Total MCMC iterations.

``--burn-in`` ``3000``
   Initial iterations discarded before sampling.

``--thin`` ``1``
   Keep one sample every ``thin`` iterations.

``--seed`` ``42``
   Random seed for reproducible MCMC.

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

   For many datasets, a practical starting point is ``--burn-in`` around
   20%-50% of ``--iters``. Increase ``--iters`` when posterior summaries are
   unstable across runs.

.. note::

   If memory is limited, reduce ``--chunk-size`` first, then enable
   ``--mmap``. This usually lowers RAM usage with a possible runtime penalty.

Examples
--------

.. code-block:: bash
   :caption: Quick Start Baseline (RR)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m RR \
      -o model_rr

Expected outputs: ``model_rr.snp.eff``, ``model_rr.param``.

.. code-block:: bash
   :caption: Accuracy-Oriented Training (R)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      -o model_bayesr

Expected outputs: ``model_bayesr.snp.eff``, ``model_bayesr.param``.

.. code-block:: bash
   :caption: Sparse Effects with Variable Selection (B)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m B \
      --pi 0.99 0.01 \
      -o model_bayesb

.. code-block:: bash
   :caption: Add Fixed Effects (qcovar + dcovar)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --dcovar sex.tsv \
      --qcovar age.tsv \
      -o model_covar

.. code-block:: bash
   :caption: Longer MCMC for Stable Posterior Estimates

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --iters 50000 \
      --burn-in 10000 \
      --thin 5 \
      -o model_high_prec

.. code-block:: bash
   :caption: Additive + Dominance Model

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      --mode AD \
      --dscale 0.0001 0.001 0.01 0.1 1.0 \
      --dpi 0.95 0.05 \
      -o model_dom

.. code-block:: bash
   :caption: Estimate Mixture Proportions

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m C \
      --pi 0.9 0.1 \
      --sample-pi \
      -o model_cpi

See Also
--------

- :ref:`predict-command` for applying trained effects to target samples.
- :ref:`assoc-command` for SNP-wise association analysis.
- :ref:`grm-command` for genomic relationship matrix construction.
