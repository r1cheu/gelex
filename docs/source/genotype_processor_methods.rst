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

If you are unsure, use ``OrthStandardizeHWE`` (default for GRM and fit). You can also use the short alias ``OSH``.

- Use ``OrthStandardizeHWE`` (alias: ``OSH``) for the default GRM/fit baseline.
- Use ``StandardizeHWE`` (alias: ``SH``) when orthogonal dominance coding is not needed.
- Use center methods (``CenterHWE``, ``OrthCenterHWE``, ``Center``, ``OrthCenter``) when only centering is needed.
- Use orthogonal methods (``OrthStandardizeHWE``, ``OrthCenterHWE``, ``OrthStandardize``, ``OrthCenter``) when orthogonal dominance coding is required.
- Use HWE methods (``StandardizeHWE``, ``CenterHWE``, ``OrthStandardizeHWE``, ``OrthCenterHWE``) when you prefer population-genetics expectations.
- Use sample methods (``Standardize``, ``Center``, ``OrthStandardize``, ``OrthCenter``) for data-driven moments.

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

- HWE methods: ``StandardizeHWE``, ``CenterHWE``, ``OrthStandardizeHWE``, ``OrthCenterHWE`` (alias: ``SH``, ``CH``, ``OSH``, ``OCH``): HWE-based expected moments
- Sample methods: ``Standardize``, ``Center``, ``OrthStandardize``, ``OrthCenter`` (alias: ``S``, ``C``, ``OS``, ``OC``): moments estimated directly from your sample

Method Matrix (User View)
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 38 32 15 15

   * - Method (alias)
     - Best for
     - Moments
     - Scaling
   * - ``StandardizeHWE`` (``SH``)
     - HWE standardization, no orthogonal dominance
     - HWE
     - Standardize
   * - ``CenterHWE`` (``CH``)
     - HWE centering, no variance scaling
     - HWE
     - Center
   * - ``OrthStandardizeHWE`` (``OSH``)
     - Default GRM/fit: orthogonal dominance + HWE
     - HWE
     - Standardize
   * - ``OrthCenterHWE`` (``OCH``)
     - Default assoc: orthogonal dominance + HWE centering
     - HWE
     - Center
   * - ``Standardize`` (``S``)
     - Sample-based standardization
     - Sample
     - Standardize
   * - ``Center`` (``C``)
     - Sample-based centering, no scaling
     - Sample
     - Center
   * - ``OrthStandardize`` (``OS``)
     - Orthogonal dominance + sample standardization
     - Sample
     - Standardize
   * - ``OrthCenter`` (``OC``)
     - Orthogonal dominance + sample centering
     - Sample
     - Center

Practical Recommendations
-------------------------

- Start with ``OrthStandardizeHWE`` (alias: ``OSH``) for most production runs.
- If biological interpretability is your top priority, prefer HWE methods
  (``StandardizeHWE``, ``CenterHWE``, ``OrthStandardizeHWE``, ``OrthCenterHWE``).
- If optimizer stability and convergence speed are your top priority, test
  sample methods (``Standardize``, ``Center``, ``OrthStandardize``, ``OrthCenter``) first.
- If comparing with older centered pipelines, use ``CenterHWE`` (alias: ``CH``).
- Use sample methods (``Standardize``, ``Center``, ``OrthStandardize``, ``OrthCenter``)
  only when you intentionally want sample-dependent centering and variance.
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

   # Recommended default (OrthStandardizeHWE)
   gelex grm -b genotypes --geno-method OSH -o grm_orth_hwe

.. code-block:: bash

   # HWE centering, no scaling
   gelex grm -b genotypes --geno-method CH -o grm_center_hwe

.. code-block:: bash

   # Orthogonal dominance with sample moments
   gelex grm -b genotypes --dom --geno-method OS -o grm_orth_sample

See Also
--------

- :ref:`grm-command`
- :doc:`api_reference`
