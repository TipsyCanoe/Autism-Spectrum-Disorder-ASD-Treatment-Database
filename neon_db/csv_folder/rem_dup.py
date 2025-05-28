import pandas as pd

# Load the CSV file
file_path = "clean_VectorExamples.csv"  # Replace with the path to your CSV file
df = pd.read_csv(file_path)

# Check if the column exists
embedding_column = "vector"  # Correct column name for embeddings
if embedding_column not in df.columns:
    raise KeyError(f"The column '{embedding_column}' does not exist in the CSV file. Available columns: {df.columns.tolist()}")

# Convert embeddings to tuples for comparison
df[embedding_column] = df[embedding_column].apply(eval)  # Convert string representation to list
df[embedding_column] = df[embedding_column].apply(tuple)  # Convert list to tuple for hashability

# Drop duplicate embeddings
df = df.drop_duplicates(subset=[embedding_column])

# Optionally, convert embeddings back to lists
df[embedding_column] = df[embedding_column].apply(list)

# Save the cleaned CSV file
output_path = "Cleaned_VectorExamples.csv"  # Replace with your desired output file name
df.to_csv(output_path, index=False)

print(f"Duplicates removed. Cleaned file saved to {output_path}.")