import click

from .convert import convert_gene_cli
from .build_lookup import build_lookup_from_fastas_cli


@click.group(invoke_without_command=True, no_args_is_help=True)
@click.version_option(version=1.0)
def entry_point():
    """Convert TCR gene names between 10X, Adaptive, and IMGT formats"""
    pass


entry_point.add_command(convert_gene_cli)
entry_point.add_command(build_lookup_from_fastas_cli)
