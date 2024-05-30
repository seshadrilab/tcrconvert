Getting started
===============

Installation
------------

``tcr-converter`` runs on Python 3 and requires ``pandas``.

You can install using ``pip``. Clone this repo, then from the top-level folder run:

.. code-block:: console

   $ pip install .

The lookup tables for translating gene names come pre-built from IMGT fasta files located under ``tcrconverter/data/``. You can re-build these tables by running:

.. code-block:: console

   $ python tcrconverter/build_lookup.py

Quick Start
-----------

Load your 10X, Adaptive, or IMGT-formatted TCR data into a ``pandas`` dataframe. It may have come from files such as:

* 10X: ``filtered_contig_annotations.csv``
* Adaptive: ``SAMPLE_TCRB.tsv``
* IMGT: A custom CSV file

Then, convert to your desired format:

.. code-block:: python

    import pandas as pd
    from tcrconverter import convert

    df = pd.read_csv("filtered_contig_annotations.csv")
    converted = convert.convert_tcr(df, fmt_from='tenx', fmt_to='adaptive')
