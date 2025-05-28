import pandas as pd
import ast

# Load the CSV file
input_file = "VectorExamples.csv"
output_file = "Cleaned_VectorExamples.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(input_file)

# Select relevant columns for the treatments table
columns_to_keep = [
    "pmid", "doi", "title", "pub_date", "abstract", "authors", "journal",
    "keywords", "url", "affiliations", "Treatment name", "Duration",
    "Primary Outcome Area", "Primary Outcome Measures", "Results: Primary measure", "vector"
]
df = df[columns_to_keep]

# Rename columns to match the database schema
df.rename(columns={
    "Treatment name": "treatment_name",
    "Duration": "duration",
    "Primary Outcome Area": "primary_outcome_area",
    "Primary Outcome Measures": "primary_outcome_measures",
    "Results: Primary measure": "results_primary_measure",
}, inplace=True)

# Handle missing values
df.fillna({
    "treatment_name": "",
    "duration": "",
    "primary_outcome_area": "",
    "primary_outcome_measures": "",
    "results_primary_measure": "",
    "vector": "[]",  # Default to an empty list for vectors
}, inplace=True)

# Ensure the vector column is properly formatted as a list of floats
def parse_vector(vector_str):
    try:
        # Convert the string representation of the vector to a Python list
        return list(map(float, ast.literal_eval(vector_str)))
    except (ValueError, SyntaxError):
        # Return an empty list if parsing fails
        return []

df["vector"] = df["vector"].apply(parse_vector)

# Save the cleaned DataFrame to a new CSV file
df.to_csv(output_file, index=False)

print(f"Cleaned data saved to {output_file}")