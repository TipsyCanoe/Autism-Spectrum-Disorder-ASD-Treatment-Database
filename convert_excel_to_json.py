import pandas as pd
from pathlib import Path

# Use: Add name of excel file in string with path, and sheet name in oldDf string
path = str(Path.cwd()) + '/Adjusted_ASD_Sheet_Info_V2.xlsx'
oldDf = pd.read_excel(path, 'Sheet1')

oldDf.to_json("Adjusted_ASD_Sheet_Info_V2.json")

