# Fetch study index from gwas catalog
from __future__ import annotations
import polars as pl
from enum import Enum
from pathlib import Path
import ftplib
import logging
from urllib.parse import urlparse, ParseResult
import io
import gzip
from typing import Callable, Sequence, TypeVar
from polars.expr.whenthen import ChainedThen
import seaborn as sns
import matplotlib.pyplot as plt

T = TypeVar("T")


class StudyIndexSchema(Enum):
    SUMSTAT_LOCATION = "SUMMARY STATS LOCATION"
    ACCESSION_ID = "STUDY ACCESSION"


class StudyIndex:
    def __init__(self, df: pl.LazyFrame):
        self.df = df
        self.schema = StudyIndexSchema

    @classmethod
    def from_catalog_studies(cls, path: str) -> StudyIndex:
        df = pl.scan_csv(path, separator="\t", quote_char="'")
        return cls(df=df)

    def find_study(self, study: str) -> StudyIndex:
        self.df = self.df.filter(pl.col(self.schema.ACCESSION_ID.value) == study)
        return self

    def find_sumstat(self, study: str) -> str:
        sumstat_path = self.find_study(study).df.select(pl.col(self.schema.SUMSTAT_LOCATION.value)).collect()
        return sumstat_path.item(0, 0)

    def get_sumstat(self, study: str, tmp_path: str) -> Sumstat:
        sumstat_uri = self.find_sumstat(study)
        sumstat = Sumstat.from_catalog_sumstat(sumstat_uri, tmp_path)
        return sumstat

    def collect(self) -> pl.DataFrame:
        return self.df.collect()


class Sumstat:
    def __init__(self, df: pl.LazyFrame):
        self.df = df
        self.columns = self.df.collect_schema().names()
        self.dtypes = self.df.collect_schema().dtypes()

    @classmethod
    def from_catalog_sumstat(cls, uri: str, path) -> Sumstat:
        protocol = uri.split("://")
        study = uri.split("/")[-1]
        match protocol:
            case ["http", _] | ["https", _]:
                parsed_uri = urlparse(uri.replace("https", "ftp").replace("http", "ftp"))
                logging.info("Reading from %s", str(uri))
                df = cache_sumstat(parsed_uri, path, study)
                return Sumstat(df=df)

            case _:
                raise ValueError("Incorrect uri to download summary statistics.")

    def manhattan(self) -> Sumstat:
        """Plot manhattan."""
        # make the order of variants correct
        self = self.drop_nullable_variants().sort().neglog_pvalue()
        data = self.df.select("base_pair_location", "neglog_pval", "chromosome").collect()
        sns.scatterplot(data=data, x="base_pair_location", y="neglog_pval", hue="chromosome")
        plt.show()
        return self

    def neglog_pvalue(self) -> Sumstat:
        """Add negative logarithm of p-value"""
        self.df = self.df.with_columns(neglog_pval=-1 * pl.col("p_value").log10())
        return self

    def drop_nullable_variants(self) -> Sumstat:
        # variant_field = "hm_variant_id" if "hm_variant_id" in self.columns else "variant_id"
        print(self.df.collect())
        self.df = self.df.filter(pl.col("rsid").is_not_null())
        return self

    def sort(self) -> Sumstat:
        if "chromosome" in self.columns:
            col_idx = self.columns.index("chromosome")
            if self.dtypes[col_idx] == pl.String():
                logging.info("Sorting by `chromosome`")

                def reduce_expression(
                    callable: Callable[[T, ChainedThen], ChainedThen],
                    initial_expression: ChainedThen,
                    collection: Sequence[T],
                ) -> pl.ChainedThen:
                    """Reduce initial expression."""
                    expression = initial_expression
                    for c in collection:
                        expression = callable(c, expression)

                    return expression

                def map_chromosomes(col: pl.Column) -> pl.Column:
                    """Map chromosomes to chromosome order."""
                    chromosomes = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]
                    order = {c: idx + 1 for idx, c in enumerate(chromosomes)}
                    expression = pl.when(col == chromosomes[0]).then(pl.lit(order[chromosomes[0]]))
                    initial_expression = pl.when(col == chromosomes[0]).then(order[chromosomes[0]])

                    def _then(c: str, e: ChainedThen) -> ChainedThen:
                        return e.when(col == c).then(pl.lit(order[c]))

                    expression = reduce_expression(
                        callable=_then,
                        initial_expression=initial_expression,  # type: ignore
                        collection=chromosomes[1:],
                    )

                    return expression

                self.df = self.df.with_columns(order=map_chromosomes(pl.col("chromosome"))).sort(
                    by=[pl.col("order"), pl.col("position")]
                )
            else:
                self.df = self.df.sort(by=[pl.col("chromosome"), pl.col("base_pair_location")])
        if "hm_chrom" in self.columns:
            logging.info("Sorting by `hm_chrom`")
            self.df = self.df.sort(by=[pl.col("hm_chrom"), pl.col("hm_pos")])

        return self


def cache_sumstat(uri: ParseResult, tmp_path: str, study: str) -> pl.LazyFrame:
    """Read and cache summary statistics in the tmp_path."""
    Path(tmp_path).mkdir(exist_ok=True, parents=True)
    sumstat_path = tmp_path + "/" + study + ".parquet"
    if Path(sumstat_path).exists():
        logging.info("Found sumstat under %s", sumstat_path)
        return pl.scan_parquet(sumstat_path)
    logging.info("Did not found sumstat in tmp, caching...")

    server = ftplib.FTP(host=uri.netloc)
    server.login()
    server.cwd(uri.path + "/harmonised")
    files = [f for f in server.nlst() if f.endswith(".h.tsv.gz")]
    assert len(files) == 1, f"Found {len(files)}, expected one file with .h.tsv.gz extension!"
    buffer = io.BytesIO()
    logging.info("Found %s", files[0])

    server.retrbinary(f"RETR {files[0]}", lambda x: buffer.write(x))
    content = gzip.decompress(buffer.getvalue())
    df = pl.read_csv(content, separator="\t", null_values=["NA"])
    buffer.close()
    df.write_parquet(sumstat_path)

    logging.info("Cache completed at %s", sumstat_path)
    return cache_sumstat(uri, tmp_path, study)


__all__ = ["StudyIndex"]
