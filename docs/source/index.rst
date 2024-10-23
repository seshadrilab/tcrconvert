.. tcr-converter documentation master file, created by
   sphinx-quickstart on Tue May 14 12:36:14 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TCRconvert
==========

Convert TCR annotations between 10X, Adaptive, and IMGT formats
----------------------------------------------------------------

The naming conventions for T-cell receptor (TCR) genes differ between sequencing 
platforms and the IMGT reference. For example, the naming of TCR alpha chain variable \
gene segment 1-2 allele 1 by the 10X and Adaptive platforms and in IMGT:

* 10X: TRAV1-2
* Adaptive: TCRAV01-02*01
* IMGT: TRAV1-2*01

TCRconvert enhances TCR dataset interoperability by providing reliable format conversion 
across 10X, Adaptive, and IMGT-formatted data. Unlike existing tools that limit 
conversions to only two formats or require custom objects, TCRconvert works directly 
with data frames TCRconvert saves researchers time and prevents errors from with manual conversions.

TCRconvert takes a Pandas DataFrame with at least one column of gene names as input. It produces a Pandas DataFrame with converted gene names as output.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage
   faq
   functions
   contributing
