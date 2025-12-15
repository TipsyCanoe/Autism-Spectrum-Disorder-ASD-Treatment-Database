import pandas as pd
from pathlib import Path
import os
from datetime import datetime

def xlsx_to_csvs(file_heads):
    today = datetime.now().strftime("%Y-%m-%d")
    combined_df = pd.DataFrame() # all data

    for file_head in file_heads:
        path = str(Path.cwd()) + '/' + file_head + '.xlsx'
        if not os.path.isfile(path):
            print(f"File not found: {path}; skipping")
            continue
        oldDf = pd.read_excel(path, 'Sheet1') # Read a file

        # Save individual CSV
        csv_name = f"{file_head}.csv"
        oldDf.to_csv(csv_name, index=False)

        # Add to combined DataFrame
        combined_df = pd.concat([combined_df, oldDf], ignore_index=True)

    # Write the combined CSV if we have data
    if not combined_df.empty:
        combined_df.to_csv(f"pubmed_combined_{today}.csv", index=False)
        print(f"Successfully converted to: pubmed_combined_{today}.csv")
    else:
        print("No data to combine")

if __name__ == '__main__':
    today = datetime.now().strftime("%Y-%m-%d")
    # Use: Add name of excel file in string with path, and sheet name in oldDf string
    file_heads = [
            f"pubmed_ASD_info_{today}",
            f"pubmed_papers_info_{today}",
            f"pubmed_treatment_info_{today}"
        ]
    xlsx_to_csvs(file_heads)

