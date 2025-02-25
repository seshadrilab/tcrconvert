# TCRconvert

> **Warning**: This project is in **beta stage**. It is under active development and may be unstable.

[![codecov](https://codecov.io/github/seshadrilab/tcrconvert/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/github/seshadrilab/tcrconvert)
[![tests](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml/badge.svg)](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/tcrconvert/badge/?version=latest)](https://tcrconvert.readthedocs.io/en/latest/?badge=latest)

**Rename T-cell receptor genes between 10X, Adaptive, and IMGT formats**

TCRconvert takes T-cell receptor (TCR) data containing V, D, J, and/or C genes from 10X, Adaptive, or other sequencing platforms and renames them from any of these formats to any other one:

* **10X**: TRAV1-2
* **Adaptive**: TCRAV01-02*01
* **IMGT**: TRAV1-2*01

TCRconvert works with human, mouse, and rhesus macaque data out-of-the-box, but users can also add their own species.

TCRconvert helps researchers unify TCR datasets by converting them to a standard naming convention. It is fast, reliable, and prevents errors from manual conversions. Unlike other tools that require custom objects, TCRconvert works directly with Pandas DataFrames and CSV/TSV files.

**You can use it two ways:**

**1. As a library**:
```python
import tcrconvert
import pandas as pd

# Convert gene names
tcr_file = tcrconvert.get_example_path('tenx.csv')
dat = pd.read_csv(tcr_file)[['barcode', 'v_gene' , 'd_gene', 'j_gene', 'c_gene', 'cdr3']]
tcrconvert.convert_gene(dat, frm='tenx', to='adaptive')

# Create a custom reference
tcrconvert.build_lookup_from_fastas('path/to/fasta/dir/', 'myspecies')
```

**2. As a command-line tool**:
```bash
$ tcrconvert convert -i tcrconvert/examples/tenx.csv -o adaptive.tsv --frm tenx --to adaptive # Convert gene names
$ tcrconvert build -i path/to/fasta_dir/ -s myspecies # Create a custom reference
```

For full documentation visit [tcrconvert.readthedocs.io](https://tcrconvert.readthedocs.io/en/latest/)

# Installation

Requirements:

* `python >=3.9`
* `pandas >=1.5.0`
* `click >=8.1.7`
* `platformdirs >=4.2.2`

TCRconvert runs on Windows, macOS, and Linux.

Install from GitHub using `pip`:

```
pip install git+https://github.com/seshadrilab/tcrconvert
```

Or clone this repo and from the top-level folder run:

```
pip install .
```

# Basic usage (library)

**Load some 10X data**


```python
import tcrconvert
import pandas as pd

tcr_file = tcrconvert.get_example_path('tenx.csv')
tcrs = pd.read_csv(tcr_file)[['barcode', 'v_gene' , 'd_gene', 'j_gene', 'c_gene', 'cdr3']]
tcrs
```





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>barcode</th>
      <th>v_gene</th>
      <th>d_gene</th>
      <th>j_gene</th>
      <th>c_gene</th>
      <th>cdr3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAACCTGAGACCACGA-1</td>
      <td>TRAV1-2</td>
      <td>TRBD1</td>
      <td>TRAJ12</td>
      <td>TRAC</td>
      <td>CAVMDSSYKLIF</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAACCTGAGACCACGA-1</td>
      <td>TRBV6-1</td>
      <td>TRBD2</td>
      <td>TRBJ2-1</td>
      <td>TRBC2</td>
      <td>CASSGLAGGYNEQFF</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAACCTGAGGCTCTTA-1</td>
      <td>TRBV6-4</td>
      <td>TRBD2</td>
      <td>TRBJ2-3</td>
      <td>TRBC2</td>
      <td>CASSGVAGGTDTQYF</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAACCTGAGGCTCTTA-1</td>
      <td>TRAV1-2</td>
      <td>TRBD1</td>
      <td>TRAJ33</td>
      <td>TRAC</td>
      <td>CAVKDSNYQLIW</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAACCTGAGTGAACGC-1</td>
      <td>TRBV2</td>
      <td>TRBD1</td>
      <td>TRBJ1-2</td>
      <td>TRBC1</td>
      <td>CASNQGLNYGYTF</td>
    </tr>
  </tbody>
</table>
</div>



**Convert gene names from the 10X format to the Adaptive format**


```python
new_tcrs = tcrconvert.convert_gene(tcrs, frm='tenx', to='adaptive')
new_tcrs
```

    Warning: Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.
    Warning: Adaptive only captures VDJ genes. Converted C genes will become NA.





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>barcode</th>
      <th>v_gene</th>
      <th>d_gene</th>
      <th>j_gene</th>
      <th>c_gene</th>
      <th>cdr3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAACCTGAGACCACGA-1</td>
      <td>TCRAV01-02*01</td>
      <td>TCRBD01-01*01</td>
      <td>TCRAJ12-01*01</td>
      <td>&lt;NA&gt;</td>
      <td>CAVMDSSYKLIF</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAACCTGAGACCACGA-1</td>
      <td>TCRBV06-01*01</td>
      <td>TCRBD02-01*01</td>
      <td>TCRBJ02-01*01</td>
      <td>&lt;NA&gt;</td>
      <td>CASSGLAGGYNEQFF</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAACCTGAGGCTCTTA-1</td>
      <td>TCRBV06-04*01</td>
      <td>TCRBD02-01*01</td>
      <td>TCRBJ02-03*01</td>
      <td>&lt;NA&gt;</td>
      <td>CASSGVAGGTDTQYF</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAACCTGAGGCTCTTA-1</td>
      <td>TCRAV01-02*01</td>
      <td>TCRBD01-01*01</td>
      <td>TCRAJ33-01*01</td>
      <td>&lt;NA&gt;</td>
      <td>CAVKDSNYQLIW</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAACCTGAGTGAACGC-1</td>
      <td>TCRBV02-01*01</td>
      <td>TCRBD01-01*01</td>
      <td>TCRBJ01-02*01</td>
      <td>&lt;NA&gt;</td>
      <td>CASNQGLNYGYTF</td>
    </tr>
  </tbody>
</table>
</div>


# Basic usage (command-line)

TCRconvert takes a `.csv` or `.tsv` file with at least one column of gene names as input. It produces a `.csv` or `.tsv` file with converted gene names as output.

**Inspect our input 10X data**


```
$ cat tcrconvert/examples/tenx.csv
```

    barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id
    AAACCTGAGACCACGA-1,TRUE,AAACCTGAGACCACGA-1_contig_1,TRUE,521,TRA,TRAV1-2,TRBD1,TRAJ12,TRAC,TRUE,TRUE,CAVMDSSYKLIF,TGTGCTGTGATGGATAGCAGCTATAAATTGATCTTC,1569,2,clonotype16,clonotype16_consensus_1
    AAACCTGAGACCACGA-1,TRUE,AAACCTGAGACCACGA-1_contig_3,TRUE,584,TRB,TRBV6-1,TRBD2,TRBJ2-1,TRBC2,TRUE,TRUE,CASSGLAGGYNEQFF,TGTGCCAGCAGTGGACTAGCGGGGGGCTACAATGAGCAGTTCTTC,5238,7,clonotype16,clonotype16_consensus_2
    AAACCTGAGGCTCTTA-1,TRUE,AAACCTGAGGCTCTTA-1_contig_1,TRUE,551,TRB,TRBV6-4,TRBD2,TRBJ2-3,TRBC2,TRUE,TRUE,CASSGVAGGTDTQYF,TGTGCCAGCAGTGGGGTAGCGGGAGGCACAGATACGCAGTATTTT,3846,4,clonotype26,clonotype26_consensus_2
    AAACCTGAGGCTCTTA-1,TRUE,AAACCTGAGGCTCTTA-1_contig_2,TRUE,518,TRA,TRAV1-2,TRBD1,TRAJ33,TRAC,TRUE,TRUE,CAVKDSNYQLIW,TGTGCTGTGAAGGATAGCAACTATCAGTTAATCTGG,2019,2,clonotype26,clonotype26_consensus_1
    AAACCTGAGTGAACGC-1,TRUE,AAACCTGAGTGAACGC-1_contig_1,TRUE,674,TRB,TRBV2,TRBD1,TRBJ1-2,TRBC1,TRUE,TRUE,CASNQGLNYGYTF,TGTGCCAGCAATCAGGGCCTTAACTATGGCTACACCTTC,3002,6,clonotype81,clonotype81_consensus_2


**Convert gene names from 10X to Adaptive**


```
$ tcrconvert convert \
    -i tcrconvert/examples/tenx.csv \
    -o tcrconvert/examples/tenx2adapt.tsv \
    --frm tenx \
    --to adaptive
```

    WARNING - Adaptive only captures VDJ genes. Converted C genes will become NA.
    WARNING - Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.



```
$ cat tcrconvert/examples/tenx2adapt.tsv
```

    barcode	is_cell	contig_id	high_confidence	length	chain	v_gene	d_gene	j_gene	c_gene	full_length	productive	cdr3	cdr3_nt	reads	umis	raw_clonotype_id	raw_consensus_id
    AAACCTGAGACCACGA-1	TRUE	AAACCTGAGACCACGA-1_contig_1	TRUE	521	TRA	TCRAV01-02*01	TCRBD01-01*01	TCRAJ12-01*01		TRUE	TRUE	CAVMDSSYKLIF	TGTGCTGTGATGGATAGCAGCTATAAATTGATCTTC	1569	2	clonotype16	clonotype16_consensus_1
    AAACCTGAGACCACGA-1	TRUE	AAACCTGAGACCACGA-1_contig_3	TRUE	584	TRB	TCRBV06-01*01	TCRBD02-01*01	TCRBJ02-01*01		TRUE	TRUE	CASSGLAGGYNEQFF	TGTGCCAGCAGTGGACTAGCGGGGGGCTACAATGAGCAGTTCTTC	5238	7	clonotype16	clonotype16_consensus_2
    AAACCTGAGGCTCTTA-1	TRUE	AAACCTGAGGCTCTTA-1_contig_1	TRUE	551	TRB	TCRBV06-04*01	TCRBD02-01*01	TCRBJ02-03*01		TRUE	TRUE	CASSGVAGGTDTQYF	TGTGCCAGCAGTGGGGTAGCGGGAGGCACAGATACGCAGTATTTT	3846	4	clonotype26	clonotype26_consensus_2
    AAACCTGAGGCTCTTA-1	TRUE	AAACCTGAGGCTCTTA-1_contig_2	TRUE	518	TRA	TCRAV01-02*01	TCRBD01-01*01	TCRAJ33-01*01		TRUE	TRUE	CAVKDSNYQLIW	TGTGCTGTGAAGGATAGCAACTATCAGTTAATCTGG	2019	2	clonotype26	clonotype26_consensus_1
    AAACCTGAGTGAACGC-1	TRUE	AAACCTGAGTGAACGC-1_contig_1	TRUE	674	TRB	TCRBV02-01*01	TCRBD01-01*01	TCRBJ01-02*01		TRUE	TRUE	CASNQGLNYGYTF	TGTGCCAGCAATCAGGGCCTTAACTATGGCTACACCTTC	3002	6	clonotype81	clonotype81_consensus_2

# Contributing

I welcome feedback! If you would like to resolve an issue or add improvements please submit a pull request.

# Issues

If you run into problems or need help running TCRconvert please file an issue on GitHub.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot `edu`

# Acknowledgments

This project was created with support from the Fred Hutchinson Cancer Center Translational Data Science Integrated Research Center (TDS IRC) through the 2024 Data Scientist Collaboration Grant.
