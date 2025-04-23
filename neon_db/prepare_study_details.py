import pandas as pd

# Load the cleaned study details and treatment ID mappings
study_details = pd.read_csv("study_details_cleaned.csv")
treatment_ids = pd.read_csv("treatment_ids.csv")

# Merge to get the treatment_id
merged = study_details.merge(treatment_ids, on="treatment_name", how="left")

# Rename 'id' column from treatment_ids.csv to 'treatment_id'
merged = merged.rename(columns={"id": "treatment_id"})

# Reorder columns so treatment_id is first
cols = ["treatment_id"] + [col for col in merged.columns if col != "treatment_id"]
final = merged[cols]

# 🔥 Drop treatment_name (not part of SQL table)
final = final.drop(columns=["treatment_name"])

# Drop rows where treatment_id is missing (no match)
final = final.dropna(subset=["treatment_id"])
final["treatment_id"] = final["treatment_id"].astype(int)

# Save final import-ready CSV
final.to_csv("study_details_final.csv", index=False, quoting=1)
print("study_details_final.csv cleaned and ready.")
