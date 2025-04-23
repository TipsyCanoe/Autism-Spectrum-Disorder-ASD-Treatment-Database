import pandas as pd
import csv
from datetime import datetime

# Load original data
df = pd.read_csv("original_data.csv")

### 1. Treatments Table ###
# Clean sample_size to keep only numbers
def clean_sample_size(val):
    try:
        return int(str(val).split()[0])  # Get the first part before space/(
    except:
        return None

treatments = df[['treatment_name', 'age_range', 'duration', 'sample_size']].copy()
treatments = treatments[treatments['treatment_name'].notnull() & treatments['duration'].notnull()]

treatments['sample_size'] = treatments['sample_size'].apply(clean_sample_size)

# Export with quoting
treatments.to_csv("treatments_cleaned.csv", index=False, quoting=csv.QUOTE_ALL)


### 2. Study Details Table ###
study_details = df[['treatment_name', 'first_author', 'doi', 'commercial_affiliation',
                    'publication_date', 'journal', 'study_type', 'gender_ratio',
                    'dose_range', 'race_ethnicity', 'notes']].copy()
study_details = study_details[study_details['treatment_name'].notnull()]
study_details.to_csv("study_details_cleaned.csv", index=False, quoting=csv.QUOTE_ALL)

### 3. Outcomes Table ###
primary = df[['treatment_name', 'primary_outcome_area', 'primary_outcome_measures', 'primary_results']].copy()
primary.columns = ['treatment_name', 'area', 'measures', 'results']
primary['outcome_type'] = 'primary'

secondary = df[['treatment_name', 'secondary_outcome_area', 'secondary_outcome_measures', 'secondary_results']].copy()
secondary.columns = ['treatment_name', 'area', 'measures', 'results']
secondary['outcome_type'] = 'secondary'

outcomes = pd.concat([primary, secondary])
outcomes = outcomes[outcomes['treatment_name'].notnull()]
outcomes.to_csv("outcomes_cleaned.csv", index=False, quoting=csv.QUOTE_ALL)

### 4. Bias Assessment Table ###
bias = df[['treatment_name', 'sequence_generation_bias', 'allocation_concealment_bias',
           'outcome_assessor_blinding_bias', 'clinician_blinding_bias',
           'attrition_bias', 'reporting_bias', 'bias_notes']].copy()
bias = bias[bias['treatment_name'].notnull()]
bias.to_csv("bias_assessment_cleaned.csv", index=False, quoting=csv.QUOTE_ALL)

print("All cleaned CSVs regenerated from original_data.csv.")
