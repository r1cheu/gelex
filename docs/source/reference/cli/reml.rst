.. _reml-command:

reml
====

Estimate variance components with REML using the Average Information (AI)
algorithm.

Use this command to fit a mixed linear model and obtain variance-component
estimates, heritability ratios, and best linear unbiased predictions (BLUPs)
from GRM and other random-effect inputs.

Basic Syntax
------------

.. code-block:: bash
   :caption: Minimum Working Command

   gelex reml -p phenotypes.tsv --grm my_grm -o reml_run

.. code-block:: bash
   :caption: Full Syntax Template

   gelex reml --pheno <pheno_file> --grm <grm_prefix...> [OPTIONS]

Required inputs are the phenotype file (``--pheno``) and at least one random
effect. A random effect can be a GRM (``--grm``), a discrete random effect
(``--drand``), a quantitative random effect (``--qrand``), or an interaction
(``--interaction``).

Options
-------

.. rubric:: Quick Start Options

``-p, --pheno`` ``required``
   Phenotype TSV file in format ``FID IID trait1 ...``.

``-o, --out`` ``gelex``
   Output prefix for result files.

.. rubric:: Random Effect Inputs

At least one random effect is required.

``--grm``
   One or more GRM prefixes (``<prefix>.bin/.id``). Each prefix contributes one
   variance component.

``--drand``
   Discrete random-effect TSV (``FID IID factor1 ...``); each factor column
   becomes a variance component via a one-hot ``ZZ^T`` kernel.

``--qrand``
   One or more quantitative random-effect matrix TSVs (``FID IID value1 ...``);
   each file forms one linear-kernel ``ZZ^T`` component.

``--interaction``
   Interaction random effect ``<a>:<b>``, the rescaled Hadamard product of two
   kernels. Each operand is a loaded effect name or a GRM prefix.

.. rubric:: Input and Covariate Options

``--pheno-col`` ``0``
   0-based trait column index after ``FID``/``IID`` (first trait = 0).

``--qcovar``
   Quantitative covariate TSV in format ``FID IID covar1 ...``.

``--dcovar``
   Discrete covariate TSV in format ``FID IID factor1 ...``.

.. rubric:: Model Configuration

``--transform`` ``none``
   Apply a rank-based Inverse-Normal Transform (INT) to the phenotype:
   ``none``, ``dint`` (Direct INT), ``iint`` (Indirect INT).

``--int-offset`` ``0.375``
   Rank-INT offset parameter ``k``.

.. rubric:: Performance Options

``--max-iter`` ``100``
   Maximum AI-REML iterations.

``--tol`` ``1e-06``
   Relative tolerance for variance-component convergence.

``-t, --threads`` ``half of available CPU cores``
   Number of CPU threads to use.

Output Files
------------

After a successful run, the following files are written:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File pattern
     - Contents
   * - ``<out>.summary``
     - Fixed-effect coefficients and variance components with estimate,
       standard error, and variance ratio.
   * - ``<out>.effects``
     - Per-individual predicted effects (BLUPs): fixed and random effect
       columns plus a ``TOTAL`` column, keyed by ``FID``/``IID``.
   * - ``<out>.log``
     - Run log with configuration, dataset summary, iteration history, and
       timing.

Examples
--------

.. code-block:: bash
   :caption: Single GRM

   gelex reml \
      -p phenotypes.tsv \
      --grm my_grm \
      -o reml_run

.. code-block:: bash
   :caption: Additive and Dominance Variance Components

   gelex reml \
      -p phenotypes.tsv \
      --grm my_grm.add my_grm.dom \
      -o reml_ad

.. code-block:: bash
   :caption: With Covariates and INT Transform

   gelex reml \
      -p phenotypes.tsv \
      --grm my_grm \
      --qcovar pcs.tsv \
      --dcovar sex.tsv \
      --transform dint \
      -o reml_covar

.. code-block:: bash
   :caption: Add a Discrete Random Effect

   gelex reml \
      -p phenotypes.tsv \
      --grm my_grm \
      --drand blocks.tsv \
      -o reml_block

See Also
--------

- :ref:`grm-command` for preparing GRM inputs.
- :ref:`assoc-command` for mixed-model GWAS sharing the same random-effect
  inputs.
