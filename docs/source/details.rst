FAQ
=====

What gene columns are required for TCRconvert?
------------------------------------------------

TCRconvert expects at least one column to be converted (V, D, J, and/or C genes). Not all gene columns are required.


What happens if my input data has missing genes?
--------------------------------------------------

Genes marked as ``NA`` in the input dataframe will also appear as ``NA`` in the output. 
If any genes are not found in the lookup table (which is based on the IMGT reference), they will be converted to ``NA``.


Are gamma-delta TCRs supported?
----------------------------------

Yes, gamma and delta chain genes for human, mouse, and rhesus macaque are supported.


How does TCRconvert handle the lack of allele information from 10X?
---------------------------------------------------------------------

Since 10X does not provide allele-level information, all genes are assigned the allele ``*01``.


How does TCRconvert handle the lack of constant (C) genes from Adaptive?
--------------------------------------------------------------------------

Adaptive does not capture C genes, so all C genes in the output data will be set to ``NA``.


What column names are expected if I'm converting data from IMGT?
------------------------------------------------------------------

IMGT does not use standard column names, so TCRconvert assumes that the 10X column names (``v_gene``, ``d_gene``, ``j_gene``, ``c_gene``) are used. 
If your data uses different names, specify them using the ``frm_cols=`` parameter.


How can I use AIRR-format files?
----------------------------------

Specify the AIRR column names you're using with ``frm_cols=`` and the gene naming convention with ``frm=``.


What if I have custom column names?
-------------------------------------

If you're using non-standard column names that do not match 10X, Adaptive, or Adaptive V2 formats, specify your column names using ``frm_cols=``.


How are gene names with "OR" or "DV" handled?
-----------------------------------------------

Gene names containing "OR" or "DV" are accounted for and can be found in the ``lookup.csv`` tables.


What happens to multiple gene names like ``TCRAV01-02/12-02``?
----------------------------------------------------------------

These gene names will be converted to ``NA`` because they are not included in the IMGT reference.


Are non-human species supported?
----------------------------------

Yes, IMGT reference FASTAs and gene tables were used to include genes for Rhesus macaque and mouse. 
Mouse genes cover both "Mouse" and "Mouse C57BL/6J" as listed in IMGT.
