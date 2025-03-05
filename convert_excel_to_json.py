import pandas as pd
from pathlib import Path

path = str(Path.cwd()) + '/Adjusted_ASD_Sheet_Info.xlsx'
oldDf = pd.read_excel(path, 'Sheet1')

oldDf.to_json("Adjusted_ASD_Sheet_Info.json")

