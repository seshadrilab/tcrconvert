# TCRconvert

> **Warning**: This project is in **beta stage**. It is under active development and may be unstable.

[![codecov](https://codecov.io/github/seshadrilab/tcrconvert/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/github/seshadrilab/tcrconvert)
[![tests](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml/badge.svg)](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/tcrconvert/badge/?version=latest)](https://tcrconvert.readthedocs.io/en/latest/?badge=latest)

**Convert human T-cell receptor (TCR) annotations between 10X, Adaptive, and IMGT formats.**

The naming conventions for T-cell receptor (TCR) genes differ between sequencing 
platforms and the IMGT reference. For example, the naming of TCR alpha chain variable 
gene segment 1-2 allele 1 by the 10X and Adaptive platforms and in IMGT:

* **10X**: TRAV1-2
* **Adaptive**: TCRAV01-02*01
* **IMGT**: TRAV1-2*01

TCRconvert enhances TCR dataset interoperability by providing reliable format conversion 
across 10X, Adaptive, and IMGT-formatted data. Unlike existing tools that limit conversions 
to only two formats or require custom objects, TCRconvert works directly with data 
frames TCRconvert saves researchers time and prevents errors from with manual conversions. 

TCRconvert takes a Pandas DataFrame with at least one column of gene names as input. It produces a Pandas DataFrame with converted gene names as output.

For full documentation, visit [tcrconvert.readthedocs.io](https://tcrconvert.readthedocs.io/en/latest/)

# Installation

`tcrconvert` requires `pandas`.

You can install from GitHub using `pip`:

```
pip install git+https://github.com/seshadrilab/tcrconvert
```

Or clone this repo and from the top-level folder run:

```
pip install .
```

The lookup tables for translating gene names come pre-built from IMGT fasta files located under ``tcrconvert/data/``

# Basic usage

**Load some 10X data**


```python
import tcrconvert
import pandas as pd

tcr_file = '/Users/emmabishop/workspace/tcrconvert/tcrconvert/examples/example_10x.csv'

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
    Warning: Adaptive only captures VDJ genes, any C genes will become NA.





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

# Contributing

I welcome feedback! If you would like to resolve an issue or add improvements please submit a pull request.

# Issues

If you run into problems or need help running TCRconvert please file an issue on GitHub.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot `edu`
