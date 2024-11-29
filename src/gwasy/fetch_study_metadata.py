import logging
import requests
from typing import Any
from enum import Enum
import re
import sys
import polars as pl


class Responses(Enum):
    OK = 200
    MovedPermanently = 301
    BadRequest = 400
    NotFound = 404


class Fields(Enum):
    Associations = "associations"
    AssociationsByStudySummary = "associationsByStudySummary"
    BackgroundEfoTraits = "backgroundEfoTraits"
    EfoTraits = "efoTraits"
    Self = "self"
    Snps = "snps"
    Study = "study"
    AccessionId = "accessionId"
    DiseaseTrait = "diseaseTrait"
    FullPvalueSet = "fullPvalueSet"


def fetch_study(api_endpoint: str, study: str) -> Any:
    """Get the study from the api_endpoint. This is based on the https://www.ebi.ac.uk/gwas/rest/docs/api#_retrieve_a_study."""
    response = requests.get(api_endpoint)
    match response.status_code:
        case Responses.OK.value:
            logging.info("Successfully connected to api endpoint %s", api_endpoint)
            study_uri: str = response.json()["_links"]["studies"]["href"]
            pattern = re.compile(r"{.*}")
            safe_uri = pattern.sub("", study_uri) + f"/{study}"
            studies_response = requests.get(safe_uri)
            match studies_response.status_code:
                case Responses.OK.value:
                    return studies_response.json()
                case _:
                    logging.error("Failed to fetch study %s from %s", study, safe_uri)
                    logging.error("Status Code: %s", studies_response.status_code)
                    sys.exit(1)
        case _:
            logging.error("Failed to connect to the api_endpoint %s", api_endpoint)
            sys.exit(1)


def parse_metadata(metadata: dict) -> pl.DataFrame:
    data = {
        Fields.Associations.name: [metadata["_links"][Fields.Associations.value]["href"]],
        Fields.AssociationsByStudySummary.name: [metadata["_links"][Fields.AssociationsByStudySummary.value]["href"]],
        Fields.BackgroundEfoTraits.name: [metadata["_links"][Fields.BackgroundEfoTraits.value]["href"]],
        Fields.EfoTraits.name: [metadata["_links"][Fields.EfoTraits.value]["href"]],
        Fields.Self.name: [metadata["_links"][Fields.Self.value]["href"]],
        Fields.Snps.name: [metadata["_links"][Fields.Snps.value]["href"]],
        Fields.Study.name: [metadata["_links"][Fields.Study.value]["href"]],
        Fields.AccessionId.name: [metadata[Fields.AccessionId.value]],
        Fields.DiseaseTrait.name: [metadata[Fields.DiseaseTrait.value]],
        Fields.FullPvalueSet.name: [metadata[Fields.FullPvalueSet.value]],
    }

    return pl.from_dict(data)
