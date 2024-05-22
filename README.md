# tcr-converter
[![codecov](https://codecov.io/gh/emjbishop/tcr-converter/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/gh/emjbishop/tcr-converter)
[![tests](https://github.com/emjbishop/tcr-converter/actions/workflows/pytest.yml/badge.svg)](https://github.com/emjbishop/tcr-converter/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/tcr-converter/badge/?version=latest)](https://tcr-converter.readthedocs.io/en/latest/?badge=latest)

**Convert human T-cell receptor (TCR) annotations between 10X, Adaptive, and IMGT formats.**

The most popular TCR sequencing platforms (Adaptive and 10X) and the widely used [IMGT](https://www.imgt.org/IMGTindex/reference.php) reference all use different naming conventions for their TCR gene annotations:

* 10X: `TRAV1-2`
* Adaptive: `TCRAV01-02*01`
* IMGT: `TRAV1-2*01`

`tcr-converter` easily converts a dataframe of TCR annotations between formats.

Check out the full documentation here: [https://tcr-converter.readthedocs.io/en/latest/](https://tcr-converter.readthedocs.io/en/latest/)

# Installation

`tcr-converter` runs on Python 3 and requires `pandas`.

You can install using `pip`. Clone this repo, then from the top-level folder run:

```
pip install -e .
```

# Usage

Load your 10X, Adaptive, or IMGT-formatted TCR data into a `pandas` dataframe. It may have come from files such as:

* 10X: `filtered_contig_annotations.csv`
* Adaptive: `SAMPLE_TCRB.tsv`
* IMGT: A custom CSV file

Then, convert to your desired format:

```
import pandas as pd
from tcrconverter import convert

df = pd.read_csv("filtered_contig_annotations.csv")
converted = convert.convert_tcr(df, fmt_from='tenx', fmt_to='adaptive')
```

# Contributing

I welcome feedback! If you would like to resolve an issue or add improvements please submit a pull request.

# Issues

If you run into problems or need help running tcr-converter please file an issue on GitHub.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot `edu`
