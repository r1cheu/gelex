.. _vi-command:

vi
==

Train SNP effect models for genomic prediction using CAVI (Coordinate Ascent
Variational Inference).

CAVI provides a deterministic, fast alternative to MCMC for BayesRR. Use this
command when you need quicker turnaround on large datasets and posterior
uncertainty quantification beyond point estimates is not required.

.. note::

   CAVI currently supports BayesRR only. For other Bayesian methods
   (BayesA/B/C/R) or full posterior diagnostics, use :ref:`mcmc-command`.

Basic Syntax
------------

.. code-block:: bash

   gelex vi -b train_data -p phenotypes.tsv -m RR -o model_vi

Options
-------

``-p, --pheno`` ``required``
   Phenotype TSV file (``FID IID trait1 ...``).

``-b, --bfile`` ``required``
   PLINK binary prefix (``.bed/.bim/.fam``).

``-m, --method`` ``RR``
   Bayesian method (CAVI supports ``RR`` only).

``--mode`` ``A``
   Genetic effect mode: ``A`` (additive), ``D`` (dominance), ``AD`` (both).

``--max-iters`` ``1000``
   Maximum CAVI iterations.

``--tol`` ``1e-6``
   Convergence tolerance (relative RSS change).

``-o, --out`` ``gelex``
   Output prefix for generated files.

See Also
--------

- :ref:`mcmc-command` for full MCMC posterior sampling.
- :ref:`predict-command` for applying trained effects to target samples.
