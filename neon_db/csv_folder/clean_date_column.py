import pandas as pd
from dateutil import parser
import sys

def clean_date(val):
    """Attempt to parse a date value into a standardized format (YYYY-MM-DD)."""
    try:
        return parser.parse(str(val), fuzzy=True).date()
    except Exception:
        return None

def clean_date_column(input_file, column_name, output_file):
    """Clean and standardize the date column in a CSV file."""
    try:
        # Load the CSV file
        df = pd.read_csv(input_file)

        # Check if the column exists
        if column_name not in df.columns:
            print(f"❌ Error: Column '{column_name}' not found in the CSV file.")
            return

        # Clean the specified column
        df[column_name] = df[column_name].apply(clean_date)

        # Save the cleaned CSV file
        df.to_csv(output_file, index=False)
        print(f"✅ Saved cleaned CSV to '{output_file}' with standardized dates in column '{column_name}'.")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    # Check for command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python clean_date_column.py <input_file> <column_name> <output_file>")
    else:
        input_file = sys.argv[1]
        column_name = sys.argv[2]
        output_file = sys.argv[3]
        clean_date_column(input_file, column_name, output_file)