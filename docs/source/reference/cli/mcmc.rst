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

If you are unsure, start with ``RR`` to establish a baseline, then try
``R`` as a stronger default for production runs.

The mixture proportions of the selection methods can either be fixed or
estimated from the data. Add ``--sample-pi`` (additive), ``--sample-dpi``
(dominance), or ``--sample-jpi`` (joint, ``CD``) to sample the proportions
instead of holding them fixed; this is more adaptive but may require longer
chains for stable estimates.

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

``--pheno-col`` ``0``
   0-based trait column index after ``FID``/``IID``; the first trait is ``0``.

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Categorical covariate TSV in format ``FID IID factor1 ...``.

.. rubric:: Model Options

``--mode`` ``A``
   Genetic effect mode: ``A`` (additive), ``D`` (dominance), or ``AD``.
   ``CD`` requires ``AD``.

``--geno-method, --gm`` ``OSH``
   Genotype coding method. Accepts codes ``SH``, ``CH``, ``OSH``, ``OCH``,
   ``S``, ``C``, ``OS``, ``OC``, ``NS``, ``NC``. See
   :ref:`genotype-processor-methods` for what each code means.

``--h2``
   Additive heritability, in the open interval ``(0, 1)``.

``--d2``
   Dominance heritability, in the open interval ``(0, 1)``.

``--random-pve``
   Variance fraction for non-SNP random effects, in ``(0, 1)``.

``--scale``
   Additive variance multipliers, used by BayesR-style models (``R``).

``--pi``
   Additive mixture proportions for the selection methods (``B``/``C``/``R``).

``--dscale``
   Dominance variance multipliers for dominance-enabled BayesR models (``R``).

``--dpi``
   Dominance mixture proportions for the selection methods (``B``/``C``/``R``).

``--jpi``
   Joint allocation proportions for ``CD``:
   both-off, additive-only, dominance-only, both-on.

``--sample-pi`` / ``--sample-dpi`` / ``--sample-jpi``
   Sample additive, dominance, or joint (``CD``) mixture proportions instead of
   holding them fixed.

``--dom-pos-prob``
   Initial probability that an active dominance effect is positive, in
   ``(0, 1)``.

.. rubric:: MCMC Options

``--iters`` ``5000``
   Total MCMC iterations.

``--burn-in`` ``3000``
   Initial iterations discarded before sampling.

``--thin`` ``1``
   Keep one sample every ``thin`` iterations.

``--seed`` ``42``
   Random seed for reproducible MCMC.

``--checkpoint-step``
   Write a checkpoint every N iterations. Omit to checkpoint only at the end.

``--from-ckpt``
   Resume the run from an existing checkpoint file.

.. rubric:: Performance

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads to use.

Output Files
------------

After a successful run, check files with your output prefix first.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - File pattern
     - Contents
     - Typical next step
   * - ``<out>.snpeff``
     - Estimated SNP effects
     - Read by ``gelex predict --gfile <out>``
   * - ``<out>.snplut``
     - Per-SNP genotype statistics used for coding
     - Read by ``gelex predict --gfile <out>``
   * - ``<out>.param``
     - Estimated fixed/covariate and random-effect terms (mean, stddev)
     - Read by ``gelex predict --gfile <out>``
   * - ``<out>.summary``
     - Posterior summary of model parameters
     - Review estimated variance components
   * - ``<out>.draws``
     - Binary posterior draws recorded during sampling
     - Read by ``gelex post --in <out>``
   * - ``<out>.log``
     - Run log and configuration used
     - Review convergence and settings

Warnings and Notes
------------------

.. note::

   For many datasets, a practical starting point is ``--burn-in`` around
   20%-50% of ``--iters``. Increase ``--iters`` when posterior summaries are
   unstable across runs.

Examples
--------

.. code-block:: bash
   :caption: Quick Start Baseline (RR)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m RR \
      -o model_rr

Expected outputs: ``model_rr.snpeff``, ``model_rr.snplut``, ``model_rr.param``, ``model_rr.summary``, ``model_rr.log``.

.. code-block:: bash
   :caption: Accuracy-Oriented Training (R)

   gelex mcmc \
      -b train_data \
      -p phenotypes.tsv \
      -m R \
      -o model_bayesr

Expected outputs: ``model_bayesr.snpeff``, ``model_bayesr.snplut``, ``model_bayesr.param``, ``model_bayesr.summary``, ``model_bayesr.log``.

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
