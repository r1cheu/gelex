.. gelex documentation master file, created by
   sphinx-quickstart on Sun Sep 21 19:04:12 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Gelex: High-Performance Genomic Analysis
========================================

.. image:: ../images/gelex_logo.jpeg
   :align: center
   :width: 100%
   :alt: Gelex Logo

*Gelex* is a high-performance C++ library and CLI tool for genomic prediction and genome-wide association studies (GWAS). It integrates advanced Bayesian models (BayesAlphabet series) and frequentist approaches (GBLUP) with memory-mapped genotype data, delivering state-of-the-art performance for large-scale genomic datasets.

.. admonition:: Quick Links
   :class: tip

   - :doc:`getting_started/installation` - Get Gelex running on your system.
   - :doc:`getting_started/quickstart` - Fit your first model in one command.
   - :doc:`tutorials/gwas` - Step-by-step guide to running your first GWAS.
   - :doc:`reference/cli/index` - Comprehensive command-line reference.

.. note::

   This project is under active development. APIs and features are subject to change.

Install the latest version via **pixi** (recommended) or **conda**, then head to
the :doc:`getting_started/quickstart`:

.. code-block:: bash

   pixi global install -c conda-forge -c https://prefix.dev/gelex gelex

Key Features
------------

*   **Bayesian Models**: Full BayesAlphabet suite (RR, A, B, C, R, CD) with dominance effect variants.
*   **Frequentist Models**: GBLUP with REML-based variance component estimation.
*   **GWAS**: Mixed linear model-based association testing with LOCO correction.
*   **High Performance**: AVX512/AVX2 vectorized I/O, OpenMP parallel processing, and optimized MKL/OpenBLAS backends.
*   **Memory Efficiency**: Memory-mapped BED file reading with chunk-based processing.

How This Documentation Is Organized
-----------------------------------

The documentation is split into four parts by what you need:

*   :doc:`Getting Started <getting_started/installation>` — install Gelex and
    run your first command.
*   **Tutorials** — learn a complete workflow end to end:
    :doc:`tutorials/genomic_selection` and :doc:`tutorials/gwas`.
*   **Concepts** — understand the statistics and design choices behind the
    tools: :doc:`Bayesian models <concepts/bayesian_models>`,
    :doc:`mixed-model association <concepts/mixed_model_gwas>`, and
    :doc:`genotype coding methods <concepts/geno_methods>`.
*   **Reference** — look up exact command options and file layouts:
    :doc:`reference/cli/index` and :doc:`reference/data_formats`.

New to Gelex? Start with :doc:`getting_started/quickstart`, then follow the
tutorial matching your goal — :doc:`tutorials/genomic_selection` for prediction
or :doc:`tutorials/gwas` for association mapping.

.. only:: not latex

   Citing Gelex
   ------------

   .. admonition:: Citation
      :class: note

      Please use the following BibTeX template to cite Gelex in scientific discourse:

      .. code-block:: bibtex

          @misc{gelex,
             author = {RuLei Chen},
             year = {2026},
             note = {https://github.com/r1cheu/gelex},
             title = {Gelex: A high-performance C++ genomic analysis toolkit}
          }

.. only:: latex

   .. rubric:: How to cite this project?

   Please use the following BibTeX template to cite Gelex in scientific discourse:

   .. code-block:: bibtex

       @misc{gelex,
          author = {RuLei Chen},
          year = {2026},
          note = {https://github.com/r1cheu/gelex},
          title = {Gelex: A high-performance C++ genomic analysis toolkit}
       }

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/genomic_selection
   tutorials/gwas

.. toctree::
   :maxdepth: 2
   :caption: Concepts

   concepts/bayesian_models
   concepts/mixed_model_gwas
   concepts/geno_methods

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/cli/index
   reference/data_formats
