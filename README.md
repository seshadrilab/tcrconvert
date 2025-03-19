# TCRconvert

[![codecov](https://codecov.io/github/seshadrilab/tcrconvert/graph/badge.svg?token=BA25XH6BS2)](https://codecov.io/github/seshadrilab/tcrconvert)
[![tests](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml/badge.svg)](https://github.com/seshadrilab/tcrconvert/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/tcrconvert/badge/?version=latest)](https://tcrconvert.readthedocs.io/en/latest/?badge=latest)

`TCRconvert` converts V, D, J, and/or C gene names between the 10X
Genomics, Adaptive Biotechnologies, and IMGT nomenclatures. It supports
alpha-beta and gamma-delta T cell receptors (TCRs) for human, mouse, and
rhesus macaque. Users can also define custom species (see the docs). An [R
version](https://github.com/seshadrilab/tcrconvertr) is also available.

Use `TCRconvert` as a Python library or on the command line.

## Background

TCR annotation tools use different gene naming conventions, making
cross-dataset searches difficult (e.g., identifying 10X-annotated TCRs
in Adaptive data). Manual conversion is complex and error-prone due to
inconsistencies in naming rules.

`TCRconvert` automates this process efficiently and accurately. Our
approach is based on analyzing multiple 10X and Adaptive data sets to
capture their naming variations.

## Installation

Requirements:

* `python >=3.9`
* `pandas >=1.5.0`
* `click >=8.1.7`
* `platformdirs >=4.2.2`

Install from GitHub:

```
pip install git+https://github.com/seshadrilab/tcrconvert
```

## Library usage

#### 1. Load TCRs into a data frame

Examples of files you may want to load:

- **10X**: `filtered_contig_annotations.csv`
- **Adaptive**: `Sample_TCRB.tsv`
- **IMGT**: Output from `MiXCR` or other tools

``` python
import tcrconvert
import pandas as pd

tcr_file = tcrconvert.get_example_path('tenx.csv')
tcrs = pd.read_csv(tcr_file)[['barcode', 'v_gene' , 'j_gene', 'cdr3']]
tcrs
#               barcode        v_gene   j_gene             cdr3
# 0  AAACCTGAGACCACGA-1     TRAV29DV5   TRAJ12     CAVMDSSYKLIF
# 1  AAACCTGAGACCACGA-1  TRBV20/OR9-2  TRBJ2-1  CASSGLAGGYNEQFF
# 2  AAACCTGAGGCTCTTA-1         TRDV2    TRDJ3  CASSGVAGGTDTQYF
# 3  AAACCTGAGGCTCTTA-1         TRGV9    TRGJ1     CAVKDSNYQLIW
```

#### 2. Convert

```python
new_tcrs = tcrconvert.convert_gene(tcrs, frm = "tenx", to = "adaptive")
#> Warning in convert_gene(tcrs, frm = "tenx", to = "adaptive"): Adaptive captures
#> only VDJ genes; C genes will be NA.
#> Converting from 10X. Using *01 as allele for all genes.
new_tcrs
#>              barcode             v_gene        j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1      TCRAV29-01*01 TCRAJ12-01*01    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TCRBV20-or09_02*01 TCRBJ02-01*01 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1      TCRDV02-01*01 TCRDJ03-01*01 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1      TCRGV09-01*01 TCRGJ01-01*01    CAVKDSNYQLIW
```

## Command-line usage

#### Use `convert` subcommand

```bash
$ tcrconvert convert --input tcrconvert/examples/tenx.csv --output adaptive.tsv --frm tenx --to adaptive
```

* `--input`: Input file path (CSV or TSV)
* `--output`: Output file path (CSV or TSV)
* `--frm`: Input TCR gene format (`tenx`, `adaptive`, `adaptivev2`, or `imgt`)
* `--to`: Output TCR gene format (`tenx`, `adaptive`, `adaptivev2`, or `imgt`)

## Contributing

Contributions are welcome! To contribute, submit a pull request. See the
[contributing page](https://tcrconvert.readthedocs.io/en/latest/contributing.html) 
for details.

## Issues

To report a bug or request a feature please open an
[issue](https://github.com/seshadrilab/tcrconvert/issues).

## Contact

For other inquiries, contact Emma Bishop: emmab5 at uw dot edu.

## Acknowledgments

This project was supported by the Fred Hutchinson Cancer Center
Translational Data Science Integrated Research Center (TDS IRC) through
the 2024 Data Scientist Collaboration Grant. Special thanks to Scott
Chamberlain for development support and Shashidhar Ravishankar for gene
name curation.
