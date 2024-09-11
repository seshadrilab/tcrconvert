# TCRconvert

[![codecov](https://codecov.io/github/seshadrilab/tcrconvert/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/github/seshadrilab/tcrconvert)
[![tests](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml/badge.svg)](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/tcrconvert/badge/?version=latest)](https://tcrconvert.readthedocs.io/en/latest/?badge=latest)

**Convert human T-cell receptor (TCR) annotations between 10X, Adaptive, and IMGT formats.**

The most popular TCR sequencing platforms (Adaptive and 10X) and the widely used [IMGT](https://www.imgt.org/IMGTindex/reference.php) reference all use different naming conventions for their TCR gene annotations. These are all the same gene:

* 10X: `TRAV1-2`
* Adaptive: `TCRAV01-02*01`
* IMGT: `TRAV1-2*01`

`tcrconvert` easily converts TCR annotations between formats.

Check out the full documentation here: [tcrconvert.readthedocs.io](https://tcrconvert.readthedocs.io/en/latest/)

# Installation

`tcrconvert` requires `pandas`.

You can install using `pip`. Clone this repo, then from the top-level folder run:

```
pip install .
```

The lookup tables for translating gene names come pre-built from IMGT fasta files located under ``tcrconvert/data/``

# Usage

Load your 10X, Adaptive, or IMGT-formatted TCR data into a `pandas` dataframe. It may have come from files such as:

* 10X: `filtered_contig_annotations.csv`
* Adaptive: `SAMPLE_TCRB.tsv`
* IMGT: A custom CSV file

Then, convert to your desired format:

```
import pandas as pd
import tcrconvert

df = pd.read_csv("filtered_contig_annotations.csv")
converted = tcrconvert.convert_gene(df, frm='tenx', to='adaptive')    # From 10X to Adaptive
```

Availble formats to convert to or from are: `tenx`, `adaptive`, `adaptivev2`, or `imgt`


To use Adaptive Immune Receptor Repertoire (AIRR) formatted files, simply specify the column names:

```
import pandas as pd
import tcrconvert

df = pd.read_csv("airr_rearrangement.tsv", sep='\t')
converted = tcrconvert.convert_gene(df, frm='tenx', to='adaptive', 
                                    frm_cols=['v_call', 'd_call', 'j_call', 'c_call'])
```

# Contributing

I welcome feedback! If you would like to resolve an issue or add improvements please submit a pull request.

# Issues

If you run into problems or need help running TCRconvert please file an issue on GitHub.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot `edu`
