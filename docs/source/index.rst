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

   - :doc:`installation` - Get Gelex running on your system.
   - :doc:`gwas_tutorial` - Step-by-step guide to running your first GWAS.
   - :doc:`cli/index` - Comprehensive command-line reference.

.. note::
   
   This project is under active development. APIs and features are subject to change.

Installation
------------

Install the latest version via **pixi** (recommended) or **conda**:

.. code-block:: bash

   # Using pixi (Global install)
   pixi global install -c conda-forge -c https://prefix.dev/gelex gelex

   # Using conda
   conda install -c conda-forge -c https://prefix.dev/gelex gelex

Quick Start
-----------

Here is how to fit a Bayesian model (BayesR) in one command:

.. code-block:: bash

   gelex fit \
     --bfile data/genotypes \
     --pheno data/phenotypes.tsv \
     --method R \
     --iters 10000 \
     --burnin 2000 \
     --o result/my_analysis

For more examples, check out the :doc:`gwas_tutorial`.

Key Features
------------

*   **Bayesian Models**: Full BayesAlphabet suite (A, B, C, R, RR) with dominance effect variants.
*   **Frequentist Models**: GBLUP with REML-based variance component estimation.
*   **GWAS**: Mixed linear model-based association testing with LOCO correction.
*   **High Performance**: AVX512/AVX2 vectorized I/O, OpenMP parallel processing, and optimized MKL/OpenBLAS backends.
*   **Memory Efficiency**: Memory-mapped BED file reading with chunk-based processing.

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

   installation
   gwas_tutorial
   gs_tutorial

.. toctree::
   :maxdepth: 2
   :caption: Reference

   data_formats
   cli/index
   api_reference