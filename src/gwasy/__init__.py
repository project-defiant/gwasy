# Perform some analysis upon gwas catalog summary statistics
from __future__ import annotations
import polars as pl


def read_summary_statistics(source: str) -> pl.DataFrame:
    """Read summary statistics"""
    return pl.scan_csv(source, separator="\t").collect()


__all__ = ["read_summary_statistics"]
