Quick Start
===========

This page gets you from a fresh install to a first result. It assumes Gelex is
already installed — see :doc:`installation` if not.

Verify the Install
------------------

Confirm the ``gelex`` binary is on your ``PATH``:

.. code-block:: bash

   gelex --help

Fit Your First Model
--------------------

Fit a Bayesian model (BayesR) on a training set:

.. code-block:: bash

   gelex mcmc \
     --bfile data/genotypes \
     --pheno data/phenotypes.tsv \
     --method R \
     --iters 10000 \
     --burn-in 2000 \
     -o result/my_analysis

This writes the fitted marker effects and summaries under the
``result/my_analysis`` prefix, ready to feed into ``gelex predict``.

Where to Go Next
----------------

- :doc:`/tutorials/genomic_selection` — full train-then-predict workflow.
- :doc:`/tutorials/gwas` — mixed-model association testing.
- :doc:`/concepts/bayesian_models` — how the BayesAlphabet methods differ.
- :doc:`/reference/cli/index` — every command and option.
- :doc:`/reference/data_formats` — input and output file layouts.
