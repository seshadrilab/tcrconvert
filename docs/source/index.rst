.. tcr-converter documentation master file, created by
   sphinx-quickstart on Tue May 14 12:36:14 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TCRconvert
==========

Convert TCR annotations between 10X, Adaptive, and IMGT formats
----------------------------------------------------------------

The most popular TCR sequencing platforms (Adaptive and 10X) and the widely used [IMGT](https://www.imgt.org/IMGTindex/reference.php) reference all use different naming conventions for their TCR gene annotations. These are all the same gene:

* 10X: `TRAV1-2`
* Adaptive: `TCRAV01-02*01`
* IMGT: `TRAV1-2*01`

`tcrconvert` easily converts TCR annotations between formats.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage
   datahandling
   contributing
   troubleshooting
   functions
