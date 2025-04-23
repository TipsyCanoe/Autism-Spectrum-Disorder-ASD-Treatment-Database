import pandas as pd

# Load cleaned outcomes and treatment ID mappings
outcomes = pd.read_csv("outcomes_cleaned.csv")
treatment_ids = pd.read_csv("treatment_ids.csv")

# Merge on treatment_name
merged = outcomes.merge(treatment_ids, on="treatment_name", how="left")

# Drop any entries without a matched ID
merged = merged.dropna(subset=["id"])

# Convert to int (removing .0)
merged["treatment_id"] = merged["id"].astype(int)

# Reorder columns
final = merged[["treatment_id", "outcome_type", "area", "measures", "results"]]

# Save final CSV
final.to_csv("outcomes_final.csv", index=False, quoting=1)  # QUOTE_ALL
print("outcomes_final.csv recreated with integer treatment IDs.")
