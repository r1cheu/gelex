.. _genotype-processor-methods:

Genotype Processing Methods
===========================

This page explains how to choose ``--geno-method`` from a user perspective.
You can use it as a quick decision guide when building GRMs.

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
- In iterative model fitting (for example, ``gelex fit``), ``*-sample``
  methods can sometimes improve numerical precision and speed up convergence.
- Use one family consistently within the same analysis workflow to avoid
  method-induced scale differences.

Quick Selection Guide
---------------------

If you are unsure, start with ``standardize``.

- Use ``standardize`` when you want the default, robust baseline.
- Use ``center`` when you want centered values without variance scaling.
- Use ``orth-*`` when your workflow explicitly requires orthogonal dominance
  coding.
- Use ``*-sample`` when you prefer sample-estimated moments over HWE-based
  moments.
- For fitting-focused workflows, consider starting with
  ``standardize-sample``.

Method Families
---------------

Scaling family:

- ``standardize*``: center and scale to unit-like variance
- ``center*``: center only

Encoding family:

- non-``orth``: dominant coding ``[0, 1, 0]``
- ``orth``: dominant coding ``[0, 2p, 4p-2]``

Orthogonal vs Non-orthogonal Dominance
--------------------------------------

- ``orth`` (``[0, 2p, 4p-2]``) is designed so additive and dominance parts are
  orthogonal under the model assumptions.
- With ``orth`` coding, whether you include dominance or not, the additive
  estimate keeps the interpretation of breeding value.
- non-``orth`` (``[0, 1, 0]``) keeps a more direct biological additive
  interpretation for the additive effect.
- Choose one coding and keep it fixed across comparable analyses.

Moment family:

- default methods: HWE-based expected moments
- ``-sample`` methods: moments estimated directly from your sample

Method Matrix (User View)
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 30 24 22

   * - CLI method
     - Best for
     - Moments source
     - Scaling
   * - ``standardize``
     - Default GRM analysis
     - HWE
     - Standardize
   * - ``center``
     - Keep original scale after centering
     - HWE
     - Center
   * - ``orth-standardize``
     - Orthogonal dominance workflows
     - HWE
     - Standardize
   * - ``orth-center``
     - Orthogonal coding without scaling
     - HWE
     - Center
   * - ``standardize-sample``
     - Data-driven moments, standardized
     - Sample
     - Standardize
   * - ``center-sample``
     - Data-driven moments, centered only
     - Sample
     - Center
   * - ``orth-standardize-sample``
     - Orthogonal + sample moments + scaling
     - Sample
     - Standardize
   * - ``orth-center-sample``
     - Orthogonal + sample moments, no scaling
     - Sample
     - Center

Practical Recommendations
-------------------------

- Start with ``standardize`` for most production runs.
- If biological interpretability is your top priority, prefer HWE methods.
- If optimizer stability and convergence speed are your top priority, test
  ``*-sample`` methods first.
- If comparing with older centered pipelines, use ``center``.
- Use ``*-sample`` only when you intentionally want sample-dependent centering
  and variance.
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

   # Recommended default
   gelex grm -b genotypes --geno-method standardize -o grm_std

.. code-block:: bash

   # Centered version
   gelex grm -b genotypes --geno-method center -o grm_center

.. code-block:: bash

   # Orthogonal coding with sample moments
   gelex grm -b genotypes --dom --geno-method orth-standardize-sample -o grm_orth

See Also
--------

- :ref:`grm-command`
- :doc:`api_reference`
