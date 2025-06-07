import pandas as pd
from pathlib import Path
import os

def xlsx_to_json(file_heads):
    for file_head in file_heads:
        path = str(Path.cwd()) + '/' + file_head + '.xlsx'
        if not os.path.isfile(path):
            print(f"File not found: {path}; skipping")
            continue
        oldDf = pd.read_excel(path, 'Sheet1')

        json_name = file_head + '.json'
        oldDf.to_json(json_name)

if __name__ == '__main__':
    # Use: Add name of excel file in string with path, and sheet name in oldDf string
    file_heads = ['pubmed_papers_info', 'pubmed_treatment_info', 'pubmed_ASD_info']
    xlsx_to_json(file_heads)

