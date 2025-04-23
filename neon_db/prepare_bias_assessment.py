import pandas as pd

# Load the cleaned bias assessment data and treatment ID mappings
bias = pd.read_csv("bias_assessment_cleaned.csv")
treatment_ids = pd.read_csv("treatment_ids.csv")

# Merge on treatment_name to get treatment_id
merged = bias.merge(treatment_ids, on="treatment_name", how="left")

# Rename 'id' column from treatment_ids to treatment_id
merged = merged.rename(columns={"id": "treatment_id"})

# Drop rows with no match
merged = merged.dropna(subset=["treatment_id"])
merged["treatment_id"] = merged["treatment_id"].astype(int)

# Reorder to match table schema
columns = [
    "treatment_id",
    "sequence_generation_bias",
    "allocation_concealment_bias",
    "outcome_assessor_blinding_bias",
    "clinician_blinding_bias",
    "attrition_bias",
    "reporting_bias",
    "bias_notes"
]
final = merged[columns]

# Save to final CSV
final.to_csv("bias_assessment_final.csv", index=False, quoting=1)  # quoting=1 = QUOTE_ALL
print("bias_assessment_final.csv created and ready to import.")
