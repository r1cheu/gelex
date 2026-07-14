.. _genotype-processor-methods:

Genotype Processing Methods
===========================

This page explains how to choose ``--geno-method`` from a user perspective.
You can use it as a quick decision guide when building GRMs or fitting models.

.. note::

   ``--geno-method`` accepts **only the short codes** shown below
   (``SH``, ``CH``, ``OSH``, ``OCH``, ``S``, ``C``, ``OS``, ``OC``,
   ``NS``, ``NC``). The long names are descriptive labels used in this
   guide; they are not valid option values.

What ``--geno-method`` Changes
------------------------------

- The way genotype values are transformed before GRM calculation.
- Whether values are only centered or also standardized.
- Whether summary statistics come from sample data or HWE assumptions.

In practice, this affects:

- **Interpretability** of additive/dominance effects
- **Numerical scale** of the GRM
- **Sensitivity** to finite-sample noise

HWE vs Sample: Practical Difference
-----------------------------------

- **HWE methods** use population-genetics expectations, so they are usually
  easier to interpret from a biological perspective.
- **Sample methods** use moments estimated from your data, so standardized
  variants are better aligned with sample-level properties (mean close to 0,
  standard deviation close to 1).
- In iterative model fitting (for example, ``gelex mcmc``), sample-moment
  methods can sometimes improve numerical precision and speed up convergence.
- Use one family consistently within the same analysis workflow to avoid
  method-induced scale differences.

Quick Selection Guide
---------------------

If you are unsure, use ``OSH`` (the default for ``grm`` and ``mcmc``).

- Use ``OSH`` (OrthStandardizeHWE) for the default GRM/model baseline.
- Use ``SH`` (StandardizeHWE) when orthogonal dominance coding is not needed.
- Use center codes (``CH``, ``OCH``, ``C``, ``OC``) when only centering is
  needed.
- Use orthogonal codes (``OSH``, ``OCH``, ``OS``, ``OC``) when orthogonal
  dominance coding is required.
- Use HWE codes (``SH``, ``CH``, ``OSH``, ``OCH``) when you prefer
  population-genetics expectations.
- Use sample codes (``S``, ``C``, ``OS``, ``OC``) for data-driven moments.
- Use NOIA codes (``NS``, ``NC``) when additive and dominance columns should be
  empirically orthogonal in structured samples.

Method Families
---------------

Scaling family:

- ``S*`` / ``*S`` (standardize): center and scale to unit-like variance
- ``C*`` / ``*C`` (center): center only

Encoding family:

- non-orthogonal codes (``SH``, ``CH``, ``S``, ``C``): dominant coding
  ``[0, 1, 0]``
- orthogonal codes (``OSH``, ``OCH``, ``OS``, ``OC``): dominant coding
  ``[0, 2p, 4p-2]``
- NOIA codes (``NS``, ``NC``): sample-frequency coding that empirically
  orthogonalizes additive and dominance columns

Orthogonal vs Non-orthogonal Dominance
--------------------------------------

- Orthogonal coding (``[0, 2p, 4p-2]``) is designed so additive and dominance
  parts are orthogonal under the model assumptions.
- With orthogonal coding, whether you include dominance or not, the additive
  estimate keeps the interpretation of breeding value.
- Non-orthogonal coding (``[0, 1, 0]``) keeps a more direct biological additive
  interpretation for the additive effect.
- Choose one coding and keep it fixed across comparable analyses.

Moment family:

- HWE codes ``SH``, ``CH``, ``OSH``, ``OCH``: HWE-based expected moments
- Sample codes ``S``, ``C``, ``OS``, ``OC``: moments estimated directly from
  your sample
- NOIA codes ``NS``, ``NC``: genotype-frequency moments estimated from your
  sample

Method Matrix (User View)
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 12 26 32 15 15

   * - Code
     - Long name
     - Best for
     - Moments
     - Scaling
   * - ``SH``
     - StandardizeHWE
     - HWE standardization, no orthogonal dominance
     - HWE
     - Standardize
   * - ``CH``
     - CenterHWE
     - HWE centering, no variance scaling
     - HWE
     - Center
   * - ``OSH``
     - OrthStandardizeHWE
     - Default GRM/model: orthogonal dominance + HWE
     - HWE
     - Standardize
   * - ``OCH``
     - OrthCenterHWE
     - Default assoc: orthogonal dominance + HWE centering
     - HWE
     - Center
   * - ``S``
     - Standardize
     - Sample-based standardization
     - Sample
     - Standardize
   * - ``C``
     - Center
     - Sample-based centering, no scaling
     - Sample
     - Center
   * - ``OS``
     - OrthStandardize
     - Orthogonal dominance + sample standardization
     - Sample
     - Standardize
   * - ``OC``
     - OrthCenter
     - Orthogonal dominance + sample centering
     - Sample
     - Center
   * - ``NS``
     - NOIAStandardize
     - Empirical additive-dominance orthogonality + standardization
     - Sample
     - Standardize
   * - ``NC``
     - NOIACenter
     - Empirical additive-dominance orthogonality + centering
     - Sample
     - Center

Practical Recommendations
-------------------------

- Start with ``OSH`` for most production runs.
- If biological interpretability is your top priority, prefer HWE codes
  (``SH``, ``CH``, ``OSH``, ``OCH``).
- If optimizer stability and convergence speed are your top priority, test
  sample codes (``S``, ``C``, ``OS``, ``OC``) first.
- If comparing with older centered pipelines, use ``CH``.
- Use NOIA codes (``NS``, ``NC``) when empirical additive-dominance
  orthogonality matters more than HWE-based moments.
- Keep method choice fixed across comparable runs to avoid scale mismatch.

Minimal Technical Notes
-----------------------

- Missing genotypes are handled automatically.
- Variants with near-zero variance are treated as monomorphic and safely
  skipped from unstable scaling.
- For frequency-based calculations, Gelex clamps estimated frequency into
  ``[0, 1]`` for numerical stability.

Example Commands
----------------

.. code-block:: bash

   # Recommended default (OSH)
   gelex grm -b genotypes --geno-method OSH -o grm_orth_hwe

.. code-block:: bash

   # HWE centering, no scaling
   gelex grm -b genotypes --geno-method CH -o grm_center_hwe

.. code-block:: bash

   # Orthogonal dominance with sample moments
   gelex grm -b genotypes --mode D --geno-method OS -o grm_orth_sample

See Also
--------

- :ref:`grm-command`
- :ref:`mcmc-command`
