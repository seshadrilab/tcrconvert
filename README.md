# tcr-converter
[![codecov](https://codecov.io/gh/emjbishop/tcr-converter/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/gh/emjbishop/tcr-converter)
[![Documentation Status](https://readthedocs.org/projects/tcr-converter/badge/?version=latest)](https://tcr-converter.readthedocs.io/en/latest/?badge=latest)

Convert human T-cell receptor (TCR) gene names between 10X, Adaptive, and IMGT formats. 

The popular 10X and Adaptive TCR sequencing platforms use different naming conventions for TCR genes, making data interoperability a challenge. The widely used IMGT reference uses yet another convention. Additional conventions exist, such as AIRR, and are currently not covered here.

An example of how the same gene is named across formats:

* 10X: `TRAV1-2`
* Adaptive: `TCRAV01-02*01`
* IMGT: `TRAV1-2*01`

Check out the full documentation here: [https://tcr-converter.readthedocs.io/en/latest/](https://tcr-converter.readthedocs.io/en/latest/)

# Usage

To convert from 10X to Adaptive:

```
import pandas as pd
from tcrconverter import convert

dat10x = pd.read_csv("example_10x.csv")
converted = convert.convert_tcr(dat10x, fmt_from='tenx', fmt_to='adaptive')
```

# Installation

Requires `pandas`. 

To install, clone the repository and run this command from the top-level folder:

```
pip install -e .
```

# Issues

If you run into problems or need help running tcr-converter please file an issue on GitHub.

# Contributing

I welcome feedback! If you would like to resolve an issue or add improvements please submit a pull request.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot `edu`
