"""
Column documentation and Excel I/O for taxa_assignment_INFO.
Output is taxa_assignment_INFO_<API>.xlsx with sheet taxa_assignment_INFO.
"""

from __future__ import annotations

import os
import pandas as pd

SHEET_DATA = "taxa_assignment_INFO"


def read_taxa_assignment_info_dataframe(xlsx_path: str) -> pd.DataFrame:
    return pd.read_excel(xlsx_path, sheet_name=SHEET_DATA, engine="openpyxl")


def write_taxa_assignment_info_xlsx(df: pd.DataFrame, xlsx_path: str) -> str:
    """Write .xlsx with sheet taxa_assignment_INFO. xlsx_path must end in .xlsx."""
    parent = os.path.dirname(xlsx_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=SHEET_DATA, index=False)
    return xlsx_path
