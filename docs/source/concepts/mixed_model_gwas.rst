.. _mixed-model-gwas:

Mixed-Model Association
=======================

This page explains the statistical model behind Gelex's GWAS: the mixed linear
model, the LOCO strategy, and the phenotype transformations. Use it to
understand the options exposed by :ref:`assoc-command`.

Genome-wide association study (GWAS) tests the statistical association between
genetic variants (SNPs) and a phenotype across the entire genome, identifying
SNPs that may point to causal genes or regulatory regions.

Mixed Linear Model (MLM)
------------------------

Gelex uses a **mixed linear model (MLM)** to account for population structure
and cryptic relatedness:

.. math::

   \mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{x}b + \mathbf{g} + \boldsymbol{\varepsilon}

where:

- :math:`\mathbf{y}` is the phenotype vector.
- :math:`\mathbf{X}\boldsymbol{\beta}` represents fixed covariate effects.
- :math:`\mathbf{x}b` is the fixed effect of the candidate SNP.
- :math:`\mathbf{g}` is the random polygenic effect, :math:`\mathbf{g} \sim \mathcal{N}(\mathbf{0},\, \sigma_g^2 \mathbf{G})`.
- :math:`\boldsymbol{\varepsilon}` is the residual error.

**Efficiency:**
Similar to tools like GCTA, EMMAX, and GEMMA, Gelex estimates the variance
components (:math:`\sigma_g^2` and :math:`\sigma_e^2`) once using a **Null Model**
(excluding the candidate SNP). These estimates are then fixed when testing each
SNP, drastically improving computational speed.

LOCO Strategy
-------------

The **Leave-One-Chromosome-Out (LOCO)** method excludes markers on the current
chromosome :math:`k` from the GRM when testing SNPs on chromosome :math:`k`.
This avoids proximal contamination and increases power.

**Implementation in Gelex:**

Gelex requires two sets of GRMs for LOCO analysis:

1. **Global GRM:** Calculated using all SNPs across the genome.
2. **Chromosome-specific GRMs:** Calculated using SNPs from each chromosome individually.

During analysis (``--loco``), Gelex loads the global GRM and dynamically
subtracts the contribution of the specific chromosome GRM to compute the
"LOCO GRM" (:math:`G_{-k}`).

.. math::
   G_{-k} = \frac{G_{whole} \cdot K_{whole} - G_k \cdot K_k}{K_{whole} - K_k}

Phenotype Transformation
------------------------

GWAS assumes approximately normally distributed residuals. For phenotypes that
deviate from normality, Gelex provides **Inverse Normal Transformation (INT)**
options:

- **Direct INT** (``dint``): Rank-based transformation applied directly to the
  raw phenotype values.
- **Indirect INT** (``iint``): Transformation applied to the residuals after
  regressing out covariate effects, preserving covariate-phenotype relationships.

Both methods use the Blom offset :math:`k = 3/8` by default:
:math:`\Phi^{-1}\!\left(\frac{r_i - k}{n - 2k + 1}\right)`, where :math:`r_i` is
the rank of individual :math:`i`.

See Also
--------

- :ref:`assoc-command` for the full option list.
- :doc:`/tutorials/gwas` for an end-to-end association workflow.
- :ref:`grm-command` for building the GRM inputs.
