.. gelex documentation master file, created by
   sphinx-quickstart on Sun Sep 21 19:04:12 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Gelex: Genome Lex
=================================================

.. image:: ../images/gelex_logo.jpeg
   :align: center
   :width: 640px
   :alt: Gelex Logo

*Gelex* is a C++ library and command-line tool for genomic prediction and
genome-wide association studies (GWAS). It implements Bayesian whole-genome
regression (the BayesAlphabet family) alongside frequentist mixed-model
approaches.

.. admonition:: Quick Links
   :class: tip

   - :doc:`getting_started/installation` - Get Gelex running on your system.
   - :doc:`getting_started/quickstart` - Fit your first model in one command.
   - :doc:`tutorials/genome_selection` - Train a model and predict breeding values.
   - :doc:`tutorials/gwas` - Step-by-step guide to running your first GWAS.
   - :doc:`reference/cli/index` - Command and option reference.

.. note::

   This project is under active development; interfaces and outputs may change
   between releases.

Packages are published to the ``gelex`` channel on `prefix.dev <https://prefix.dev>`_.
Install via pixi or conda, then continue with the
:doc:`getting_started/quickstart`:

.. code-block:: bash

   pixi global install -c conda-forge -c https://prefix.dev/gelex gelex

Scope
-----

*   **Bayesian whole-genome regression** — the BayesAlphabet family for
    genome selection, under additive, dominance, or joint additive-dominance
    effect modes (:doc:`concepts/bayesian_models`).
*   **Variance component estimation** — GBLUP with AI-REML
    (:doc:`reference/cli/reml`).
*   **Association testing** — mixed linear models with LOCO correction
    (:doc:`concepts/mixed_model_gwas`).

How This Documentation Is Organized
-----------------------------------

The documentation is split into four parts by what you need:

*   :doc:`Getting Started <getting_started/installation>` — install Gelex and
    run your first command.
*   **Tutorials** — learn a complete workflow end to end:
    :doc:`tutorials/genome_selection` and :doc:`tutorials/gwas`.
*   **Concepts** — understand the statistics and design choices behind the
    tools: :doc:`Bayesian models <concepts/bayesian_models>`,
    :doc:`mixed-model association <concepts/mixed_model_gwas>`, and
    :doc:`genotype coding methods <concepts/geno_methods>`.
*   **Reference** — look up exact command options and file layouts:
    :doc:`reference/cli/index` and :doc:`reference/file_formats`.

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
             title = {Gelex: A C++ toolkit for genomic prediction and association studies}
          }

.. only:: latex

   .. rubric:: How to cite this project?

   Please use the following BibTeX template to cite Gelex in scientific discourse:

   .. code-block:: bibtex

       @misc{gelex,
          author = {RuLei Chen},
          year = {2026},
          note = {https://github.com/r1cheu/gelex},
          title = {Gelex: A C++ toolkit for genomic prediction and association studies}
       }

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/genome_selection
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
   reference/file_formats
