import os
import pytest
import pandas as pd
from convert_excels_to_csvs import xlsx_to_csvs

@pytest.mark.parametrize("file_head", [
    "pubmed_papers_info",
    "pubmed_treatment_info",
    "pubmed_ASD_info"
])

# Tests no csv file for populating database
def test_csv_file_no_exist(tmp_path, monkeypatch, file_head):
    monkeypatch.chdir(tmp_path)
    xlsx_to_csvs(file_head)
    for file_head in ['pubmed_papers_info', 'pubmed_treatment_info', 'pubmed_ASD_info']:
        assert not (tmp_path / f"{file_head}.csv").exists()

# Tests file exists after conversion
def test_csv_file_convert(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    file_heads = ['pubmed_papers_info', 'pubmed_treatment_info', 'pubmed_ASD_info']

    for file_head in file_heads:
        # Create dummy Excel file
        df = pd.DataFrame({"col1": [1, 2], "col2": ["a", "b"]})
        excel_path = tmp_path / f"{file_head}.xlsx"
        df.to_excel(excel_path, sheet_name="Sheet1", index=False)
        
    # Run conversion
    xlsx_to_csvs(file_heads)

    csv_path = tmp_path / f"{file_head}.csv"
    assert csv_path.exists(), f"{csv_path} should exist after conversion"