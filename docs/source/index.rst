.. tcr-converter documentation master file, created by
   sphinx-quickstart on Tue May 14 12:36:14 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TCRconvert
===========

Rename T-cell receptor genes between 10X, Adaptive, and IMGT formats
---------------------------------------------------------------------

TCRconvert takes T-cell receptor (TCR) data containing V, D, J, and/or C genes 
from 10X, Adaptive, or other sequencing platforms and renames them from any of 
these formats to any other one:

* **10X**: TRAV1-2
* **Adaptive**: TCRAV01-02*01
* **IMGT**: TRAV1-2*01

TCRconvert works with human, mouse, and rhesus macaque data out-of-the-box, but 
users can also add their own species (see "Using a custom reference" in the usage pages).

TCRconvert helps researchers unify TCR datasets by converting them to a standard 
naming convention. It is fast, reliable, and prevents errors from manual conversions. 
Unlike other tools that require custom objects, TCRconvert works directly with 
Pandas DataFrames and CSV/TSV files.

**You can use it two ways:**

**1. As a library**:

.. code-block:: python

   import tcrconvert

   tcrconvert.convert_gene(dat, frm='tenx', to='adaptive')  # Convert gene names
   tcrconvert.build_lookup_from_fastas('path/to/fasta_dir/', 'rabbit')  # Create a custom reference


**2. As a command-line tool**:

.. code-block:: console

   $ tcrconvert convert -i 10x.csv -o adaptive.tsv --frm tenx --to adaptive # Convert gene names
   $ tcrconvert build -i path/to/fasta_dir/ -s rabbit # Create a custom reference


View on `GitHub <https://github.com/seshadrilab/tcrconvert>`_.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage_library
   usage_cli
   faq
   functions
   contributing
