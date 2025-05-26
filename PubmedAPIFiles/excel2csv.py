import pandas as pd
import os

def main():
    # Prompt for Excel file
    excel_file = input("Enter the Excel filename (with or without .xlsx): ").strip()

    # Add extension if missing
    if not excel_file.lower().endswith('.xlsx'):
        excel_file += '.xlsx'

    # Check if file exists
    if not os.path.exists(excel_file):
        print(f"File not found: {excel_file}")
        return

    # Set output CSV name
    csv_file = os.path.splitext(excel_file)[0] + '.csv'

    try:
        # Read and convert
        df = pd.read_excel(excel_file)
        df.to_csv(csv_file, index=False)
        print(f"Successfully converted to: {csv_file}")
    except Exception as e:
        print(f"Error during conversion: {e}")

if __name__ == "__main__":
    main()
