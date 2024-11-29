import logging
import typer
from gwasy.fetch_study_index import StudyIndex

logging.basicConfig(level=logging.DEBUG)


def _main(
    study: str, path: str = "gwas-catalog-v1.0.3.1-studies-r2024-11-20.tsv", tmp_path: str = "/tmp/catalog/sumstat"
) -> None:
    sumstats = StudyIndex.from_catalog_studies(path).get_sumstat(study, tmp_path)
    print(sumstats.df.collect())


def main():
    typer.run(_main)
