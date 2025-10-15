import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser(description="Convert a JSON file to a CSV file using pandas.")
    parser.add_argument("basename", help="Base name of the file (without extension). Example: 'data' converts data.json -> data.csv")
    args = parser.parse_args()

    json_file = f"{args.basename}.json"
    csv_file = f"{args.basename}.csv"

    if not os.path.exists(json_file):
        print(f"Error: {json_file} not found.")
        return

    try:
        df = pd.read_json(json_file)

        df.to_csv(csv_file, index=False)
        print(f"Converted {json_file} -> {csv_file}")
    except Exception as e:
        print(f"Failed to convert {json_file}: {e}")

if __name__ == "__main__":
    main()
