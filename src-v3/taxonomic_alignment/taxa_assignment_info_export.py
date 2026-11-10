"""
Column documentation and Excel I/O for taxa_mapping_INFO.
Output is taxa_mapping_INFO_<API>.xlsx with sheet taxa_mapping_INFO.
"""

from __future__ import annotations

import os
import pandas as pd

SHEET_DATA = "taxa_mapping_INFO"


def read_taxa_assignment_info_dataframe(xlsx_path: str, **read_excel_kwargs) -> pd.DataFrame:
    # Callers that only need to edit the sheet keep the old defaults.
    # The final file check passes keep_default_na=False so text such as "NA"
    # is not treated as a missing value and then written back blank.
    options = {"sheet_name": SHEET_DATA, "engine": "openpyxl"}
    options.update(read_excel_kwargs)
    return pd.read_excel(xlsx_path, **options)


def write_taxa_assignment_info_xlsx(df: pd.DataFrame, xlsx_path: str) -> str:
    """Write .xlsx with sheet taxa_mapping_INFO. xlsx_path must end in .xlsx."""
    parent = os.path.dirname(xlsx_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=SHEET_DATA, index=False)
    return xlsx_path
