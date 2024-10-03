Details
=======

Here are some extra notes on how TCRconvert works.

Missing V, D, J, and/or C columns
---------------------------------

TCRconvert converts data for the gene columns that are in your data, not all gene columns are required.


Missing genes
-------------

Genes that are ``NA`` in the input dataframe will be ``NA`` in the output dataframe. 
Genes that are not in the lookup table (which is based on the IMGT reference) will be converted to ``NA``.


Converting from 10X
-------------------

Becuase 10X does not provide allele-level information, all genes are assigned allele ``*01``.


Converting from IMGT
--------------------

Becuase IMGT does not have standard column names, TCRconvert assumes 10X column names are used: 'v_gene', 'd_gene', 'j_gene', 'c_gene'. 
If not, please specify using ``frm_cols=``.


Converting to Adaptive
----------------------

Adaptive does not capture C genes, so all C genes become NA in the output data.


AIRR data
---------

The AIRR community specifies column names, but not gene naming convention. To use AIRR-format files, specify their column names using ``frm_cols=``.


Custom column names
-------------------

If you're not using standard 10X, Adaptive, or Adaptive V2 column names, please specify your column names using e.g. ``frm_cols=['myV', 'myD', 'myJ']``.


Gene names with "OR" or "DV"
----------------------------

These are accounted for and can be seen in the ``lookup.csv`` tables.


Multiple gene names, e.g. ``TCRAV01-02/12-02``
----------------------------------------------

These will become ``NA`` because they are not in the IMGT reference.
