import pandas as pd
from pathlib import Path
import os

def xlsx_to_csvs(file_heads):
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
        combined_df.to_csv("pubmed_combined.csv", index=False)
        print("Successfully converted to: pubmed_combined.csv")
    else:
        print("No data to combine")

if __name__ == '__main__':
    # Use: Add name of excel file in string with path, and sheet name in oldDf string
    file_heads = ['pubmed_papers_info', 'pubmed_treatment_info', 'pubmed_ASD_info']
    xlsx_to_csvs(file_heads)

