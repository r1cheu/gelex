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

If you are unsure, use ``3`` (``orth-standardize-hwe``, default for GRM and fit).

- Use ``3`` (orth-standardize-hwe) for the default GRM/fit baseline.
- Use ``1`` (standardize-hwe) when orthogonal dominance coding is not needed.
- Use even-numbered methods (``2``, ``4``, ``6``, ``8``) when only centering is needed.
- Use ``orth`` methods (``3``, ``4``, ``7``, ``8``) when orthogonal dominance coding is required.
- Use ``*-hwe`` methods (``1``–``4``) when you prefer population-genetics expectations.
- Use sample methods (``5``–``8``) for data-driven moments.

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

- Methods ``1``–``4`` (``*-hwe``): HWE-based expected moments
- Methods ``5``–``8``: moments estimated directly from your sample

Method Matrix (User View)
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 8 30 32 15 15

   * - Int
     - CLI name
     - Best for
     - Moments
     - Scaling
   * - ``1``
     - ``standardize-hwe``
     - HWE standardization, no orthogonal dominance
     - HWE
     - Standardize
   * - ``2``
     - ``center-hwe``
     - HWE centering, no variance scaling
     - HWE
     - Center
   * - ``3``
     - ``orth-standardize-hwe``
     - Default GRM/fit: orthogonal dominance + HWE
     - HWE
     - Standardize
   * - ``4``
     - ``orth-center-hwe``
     - Default assoc: orthogonal dominance + HWE centering
     - HWE
     - Center
   * - ``5``
     - ``standardize``
     - Sample-based standardization
     - Sample
     - Standardize
   * - ``6``
     - ``center``
     - Sample-based centering, no scaling
     - Sample
     - Center
   * - ``7``
     - ``orth-standardize``
     - Orthogonal dominance + sample standardization
     - Sample
     - Standardize
   * - ``8``
     - ``orth-center``
     - Orthogonal dominance + sample centering
     - Sample
     - Center

Practical Recommendations
-------------------------

- Start with ``3`` (orth-standardize-hwe) for most production runs.
- If biological interpretability is your top priority, prefer HWE methods (``1``–``4``).
- If optimizer stability and convergence speed are your top priority, test
  sample methods (``5``–``8``) first.
- If comparing with older centered pipelines, use ``2`` (center-hwe).
- Use sample methods (``5``–``8``) only when you intentionally want
  sample-dependent centering and variance.
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

   # Recommended default (orth-standardize-hwe)
   gelex grm -b genotypes --geno-method 3 -o grm_orth_hwe

.. code-block:: bash

   # HWE centering, no scaling
   gelex grm -b genotypes --geno-method 2 -o grm_center_hwe

.. code-block:: bash

   # Orthogonal dominance with sample moments
   gelex grm -b genotypes --dom --geno-method 7 -o grm_orth_sample

See Also
--------

- :ref:`grm-command`
- :doc:`api_reference`
