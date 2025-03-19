FAQ
=====

How does TCRconvert work?
---------------------------

TCRconvert essentially performs a ``merge`` between the input data and a 
lookup table that includes the naming conventions for each gene. 
These lookup tables are constructed from IMGT reference FASTA files and account 
for the specific naming peculiarities of each platform. The built-in lookup tables 
are located under ``tcrconvert/data/``.


What input columns are required?
----------------------------------

TCRconvert expects at least one V, D, J, and/or C gene column. You can use standard 
10X and Adaptive column names or custom names.


What if I have missing genes?
-------------------------------

``NA`` values in the input dataframe will remain ``NA`` in the output. Genes that 
are not found in the lookup table (which is based on the IMGT reference), will be 
converted to ``NA``. The built-in lookup tables are located under ``tcrconvert/data/``.


Are gamma-delta TCRs supported?
----------------------------------

Yes, for human, mouse, and rhesus macaque.


How are alleles added from 10X data?
--------------------------------------

Since 10X does not provide allele-level information, all genes are assigned the allele ``*01``.


How are C genes converted to Adaptive?
----------------------------------------

Adaptive does not capture constant ("C") gene information. If converting to the 
Adaptive format, all C genes will be set to ``NA``.


What column names should I use for my IMGT-formatted data?
------------------------------------------------------------

IMGT does not have standard column names, so it's assumed that the 10X names 
are used: (``v_gene``, ``d_gene``, ``j_gene``, ``c_gene``). To use other names, 
specify them as a list with ``frm_cols`` (library) or ``--column`` (command-line).


Can I input AIRR files?
-------------------------

Yes, just specify the AIRR column names (``v_call``, ``d_call``, ``j_call``, ``c_call``) 
using ``frm_cols`` (library) or ``--column`` (command-line). You must still 
specify the input naming convention with ``frm``.


What if I have custom column names?
-------------------------------------

If you're using non-standard column names that do not match 10X, Adaptive, or 
Adaptive V2 formats, specify them with ``frm_cols`` (library) or 
``--column`` (command-line).


What about odd names (e.g. ``TRAV14DV4``, ``TCRAV01-02/12-02``)?
------------------------------------------------------------------

Gene names containing "OR" or "DV" are accounted for in the lookup tables.

Combinations of gene names, like ``TCRAV01-02/12-02``, will be converted to ``NA`` 
because they are not in the IMGT reference.


Are non-human species supported?
----------------------------------

Mouse and rhesus macaque are supported out-of-the-box. For other species, see 
pages for building custom lookup tables for library or command-line usage.

The rhesus and mouse lookup tables were built from IMGT reference FASTAs and 
gene tables. Mouse genes cover both "Mouse" and "Mouse C57BL/6J" as listed in IMGT.


What if my Adaptive data lacks ``x_resolved`` / ``xMaxResolved`` columns?
---------------------------------------------------------------------------

Create them by combining ``x_gene``/``xGeneName`` and 
``x_allele``/``xGeneAllele`` with ``*`` as a separator. Example code:

.. code-block:: python


    import pandas as pd
    import numpy as np

    def create_col(gene_col, allele_col):
        return np.where(allele_col.notna(), gene_col + "*" + allele_col, gene_col)

    new_df = adaptive_df.copy()

    # Adaptive
    new_df['v_resolved'] = create_col(new_df['v_gene'], new_df['v_allele'])
    new_df['d_resolved'] = create_col(new_df['d_gene'], new_df['d_allele'])
    new_df['j_resolved'] = create_col(new_df['j_gene'], new_df['j_allele'])

    # Adaptive v2
    new_df['vMaxResolved'] = create_col(new_df['vGeneName'], new_df['vGeneAllele'])
    new_df['dMaxResolved'] = create_col(new_df['dGeneName'], new_df['dGeneAllele'])
    new_df['jMaxResolved'] = create_col(new_df['jGeneName'], new_df['jGeneAllele'])
