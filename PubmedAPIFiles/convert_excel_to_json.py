import pandas as pd
from pathlib import Path

# Use: Add name of excel file in string with path, and sheet name in oldDf string
file_heads = ['pubmed_papers_info.xlsx', 'pubmed_treatment_info.xlsx', 'pubmed_ASD_info.xlsx']

for file_head in file_heads:
    path = str(Path.cwd()) + '/' + file_head
    oldDf = pd.read_excel(path, 'Sheet1')

    json_name = file_head + '.json'
    oldDf.to_json(json_name)

