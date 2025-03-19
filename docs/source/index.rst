.. tcr-converter documentation master file, created by
   sphinx-quickstart on Tue May 14 12:36:14 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TCRconvert
============

Convert TCR gene names
------------------------

``TCRconvert`` converts V, D, J, and/or C gene names between the 10X
Genomics, Adaptive Biotechnologies, and IMGT nomenclatures. It supports
alpha-beta and gamma-delta T cell receptors (TCRs) for human, mouse, and
rhesus macaque. Users can also define custom species (see the docs). An 
`R version <https://github.com/seshadrilab/tcrconvertr>`_ is also available.

Use ``TCRconvert`` as a Python library or on the command line.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage_library
   custom_species_lib
   usage_cli
   custom_species_cli
   faq
   functions
   contributing


Background
------------

TCR annotation tools use different gene naming conventions, making
cross-dataset searches difficult (e.g., identifying 10X-annotated TCRs
in Adaptive data). Manual conversion is complex and error-prone due to
inconsistencies in naming rules.

``TCRconvert`` automates this process efficiently and accurately. Our
approach is based on analyzing multiple 10X and Adaptive data sets to
capture their naming variations.


Issues
--------

To report a bug or request a feature please open an
`issue <https://github.com/seshadrilab/tcrconvert/issues>`_.

Contact
---------

For other inquiries, contact Emma Bishop: emmab5 at uw dot edu.

Acknowledgments
-----------------

This project was supported by the Fred Hutchinson Cancer Center
Translational Data Science Integrated Research Center (TDS IRC) through
the 2024 Data Scientist Collaboration Grant. Special thanks to Scott
Chamberlain for development support and Shashidhar Ravishankar for gene
name curation.
