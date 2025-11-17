.. gelex documentation master file, created by
   sphinx-quickstart on Sun Sep 21 19:04:12 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
==================================

*Gelex* is a high-performance C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides efficient computation for genomic selection with memory-mapped genotype data support, delivering superior performance for large-scale genomic datasets.

Gelex implements state-of-the-art genomic prediction methods with a focus on computational efficiency and scalability. The library leverages modern C++23 standards, multi-threaded computation with OpenMP, and optimized linear algebra backends (MKL or OpenBLAS).

**Key Features:**

- **Bayesian Models**: Full classic BayesAlphabet suite including BayesA, BayesB(π), BayesC(π), BayesR, BayesRR with dominance effects.
- **High Performance**: utilizes OpenMP for parallel processing and optimized BLAS libraries for linear algebra operations
- **Memory Efficiency**: Memory-mapped and chunk-based processing
- **Modern Architecture**: Template-based MCMC implementation with hierarchical effect management
- **Comprehensive CLI**: Easy-to-use command-line interface for both fitting and prediction tasks

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite Gelex in scientific discourse:

.. code-block:: bibtex

    @misc{gelex,
       author = {Rui Chen},
       year = {2025},
       note = {https://github.com/r1cheu/gelexy},
       title = {Gelex: High-performance genomic prediction with Bayesian and frequentist approaches}
    }

.. only:: not latex

   Table of contents
   -----------------
.. toctree::
   :maxdepth: 2

   installing

.. toctree::
   :maxdepth: 2
   :caption: Basic

   cli_usage
   data_formats



.. toctree::
   :maxdepth: 2
   :caption: Model Reference

   bayesian_models
   frequentist_models

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   api_reference
