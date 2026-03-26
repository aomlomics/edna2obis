"""
Column documentation and Excel I/O for taxa_assignment_INFO.
Output is taxa_assignment_INFO_<API>.xlsx with sheets taxa_assignment_INFO and README.
"""

from __future__ import annotations

import os
import pandas as pd

SHEET_DATA = "taxa_assignment_INFO"
SHEET_README = "README"

# (column name, short description). WoRMS-only or GBIF-only called out in the text.
TAXA_ASSIGNMENT_COLUMN_ROWS: list[tuple[str, str]] = [
    ("verbatimIdentification", "Original taxonomic string from your data before matching."),
    ("cleanedTaxonomy", "Normalized / cleaned version of verbatimIdentification (what is actually used for API lookup."),
    ("ambiguous", "WoRMS. True when more than one candidate assignment was found for this verbatim string."),
    ("name_change", "WoRMS. True when WoRMS replaced an unaccepted name with the accepted valid name."),
    (
        "unaccepted_match_row",
        "WoRMS. True when the assignment in that row is an unaccepted name. Shows you the unaccepted name for comparison.",
    ),
    ("ranks_matched", "WoRMS. Number of taxonomic ranks given in your verbatimIdentification that perfectly matched with the WoRMS assignment's ranks."),
    ("ranks_provided", "WoRMS. Number of taxonomic ranks given in your verbatimIdentification for that taxonomy."),
    ("assignment_score", "WoRMS. Ratio of ranks_matched / ranks_provided."),
    ("environment", "WoRMS habitat labels when available. Non-marine taxa can be included in results when worms_return_all_matches is true in config.yaml."),
    ("selected_match", "True on the row chosen for Occurrence Core. For taxassign results, this is True for the match that WOULD be chosen if you ran the full pipeline."),
    ("scientificName", "Scientific name for this candidate row from the backbone."),
    ("confidence", "GBIF only. Assignment confidence score from 0 to 100."),
    ("taxonRank", "Rank of scientificName assigned."),
    ("scientificNameID", "WoRMS only. LSID for the taxon."),
    ("taxonID", "GBIF only. Backbone taxonomic identifier."),
    (
        "higherClassification",
        "WoRMS only. Pipe-separated lineage when worms_return_higher_classification is enabled in config.yaml. Slows performance.",
    ),
    ("kingdom", "Darwin Core kingdom."),
    ("phylum", "Darwin Core phylum."),
    ("class", "Darwin Core class."),
    ("order", "Darwin Core order."),
    ("family", "Darwin Core family."),
    ("genus", "Darwin Core genus."),
    ("match_type_debug", "Short code for how this row was matched via the API."),
    ("nameAccordingTo", "Which API you used to assign taxonomy."),
]


def taxa_assignment_readme_dataframe() -> pd.DataFrame:
    return pd.DataFrame(TAXA_ASSIGNMENT_COLUMN_ROWS, columns=["column", "description"])


def read_taxa_assignment_info_dataframe(xlsx_path: str) -> pd.DataFrame:
    return pd.read_excel(xlsx_path, sheet_name=SHEET_DATA, engine="openpyxl")


def write_taxa_assignment_info_xlsx(df: pd.DataFrame, xlsx_path: str) -> str:
    """Write .xlsx with sheets taxa_assignment_INFO and README. xlsx_path must end in .xlsx."""
    parent = os.path.dirname(xlsx_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    readme = taxa_assignment_readme_dataframe()
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=SHEET_DATA, index=False)
        readme.to_excel(writer, sheet_name=SHEET_README, index=False)
    return xlsx_path
