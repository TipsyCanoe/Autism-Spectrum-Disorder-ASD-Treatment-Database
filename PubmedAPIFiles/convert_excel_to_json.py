import pandas as pd
from pathlib import Path

# Use: Add name of excel file in string with path, and sheet name in oldDf string
file_heads = ['pubmed_papers_info', 'Adjusted_ASD_Sheet_Info_V2']

for file_head in file_heads:
    path = str(Path.cwd()) + '/' + file_head + '.xlsx'
    oldDf = pd.read_excel(path, 'Sheet1')

    json_name = file_head + '.json'
    oldDf.to_json(json_name)

