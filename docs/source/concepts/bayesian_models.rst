.. _bayesian-models:

Bayesian Models (BayesAlphabet)
===============================

This page explains the Bayesian models Gelex uses for genomic prediction and
how their priors differ. Use it to choose a ``--method`` for :ref:`mcmc-command`.

Why Bayesian Whole-Genome Regression
------------------------------------

Genomic Selection (GS) uses genome-wide markers to predict complex traits.
Unlike GWAS, which focuses on identifying individual significant variants, GS
aims to capture the *total* genetic value of an individual by simultaneously
estimating the effects of all markers, even those with small effects.

Gelex implements the **BayesAlphabet** family of models. They share the same
whole-genome regression backbone and differ only in the prior placed on marker
effects, which encodes an assumption about the trait's genetic architecture —
from oligogenic (few large effects) to highly polygenic (many small effects).

Choosing a Method
-----------------

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - ``--method``
     - Model
     - Prior assumption
   * - ``RR``
     - BayesRR / Ridge
     - All SNPs have non-zero effects drawn from a single normal distribution.
       Good for highly polygenic traits.
   * - ``A``
     - BayesA
     - All SNPs have non-zero effects, but each has its own variance
       (heavy-tailed), allowing some larger effects.
   * - ``B``
     - BayesB
     - Variable selection: a proportion ``pi`` of SNPs have zero effect; the
       rest follow a BayesA-like heavy-tailed prior.
   * - ``C``
     - BayesC
     - Variable selection: a proportion ``pi`` of SNPs have zero effect; the
       rest share a common effect variance.
   * - ``R``
     - BayesR
     - Flexible mixture of normal distributions with different variances
       (e.g. zero, small, medium, large). Often the highest accuracy across
       diverse architectures, and a good default.
   * - ``CD``
     - BayesCD
     - Couples additive and dominance marker allocation in a joint mixture.
       Requires ``--mode AD``; more parameters and longer runtime.

.. tip::
   If you are unsure, start with **BayesR** (``--method R``). It adapts to a wide
   range of genetic architectures without committing to a single assumption.

See Also
--------

- :ref:`mcmc-command` for the full option list.
- :doc:`/tutorials/genomic_selection` for an end-to-end prediction workflow.
