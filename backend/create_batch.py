import pandas as pd
import json
from openai import OpenAI

client = OpenAI()

df = pd.read_csv('grabbed_papers.csv')
batch_size = 5

schema = {
    "type": "object",
    "properties": {
        "results": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "pmid": {"type": "string"},
                    "sample_size": {"type": ["number", "null"]},
                    "treatment_name": {"type": ["string", "null"]},
                    "study_type": {"type": ["string", "null"]},
                    "study_duration": {"type": ["string", "null"]},
                    "male_female_ratio": {"type": ["string", "null"]},
                    "age_range": {"type": ["string", "null"]},
                    "treatment_dose_range": {"type": ["string", "null"]},
                    "primary_outcome_area": {"type": ["string", "null"]},
                    "primary_outcome_measures": {"type": ["string", "null"]},
                    "results_primary_measure": {"type": ["string", "null"]},
                    "secondary_outcome_area": {"type": ["string", "null"]},
                    "secondary_outcome_measures": {"type": ["string", "null"]},
                    "results_secondary_measures": {"type": ["string", "null"]},
                    "tolerability_side_effects": {"type": ["string", "null"]},
                    "safety": {"type": ["string", "null"]},
                    "drop_out_rate": {"type": ["string", "null"]},
                    "ethnicity_percentages": {"type": ["string", "null"]}
                },
                "required": [
                    "pmid",
                    "sample_size",
                    "treatment_name",
                    "study_type",
                    "study_duration",
                    "male_female_ratio",
                    "age_range",
                    "treatment_dose_range",
                    "primary_outcome_area",
                    "primary_outcome_measures",
                    "results_primary_measure",
                    "secondary_outcome_area",
                    "secondary_outcome_measures",
                    "results_secondary_measures",
                    "tolerability_side_effects",
                    "safety",
                    "drop_out_rate",
                    "ethnicity_percentages"
                ],
                "additionalProperties": False
            }
        }
    },
    "required": ["results"],
    "additionalProperties": False
}

with open('batch_requests.jsonl', 'w', encoding='utf-8') as f:
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size]
        
        papers = []
        for _, row in batch.iterrows():
            papers.append({
                "pmid": str(row['pmid']) if pd.notna(row['pmid']) else 'N/A',
                "abstract": row['abstract'] if pd.notna(row['abstract']) else 'N/A'
            })
        
        json_object = {
            "custom_id": f"batch_{i//batch_size}",
            "method": "POST",
            "url": "/v1/chat/completions",
            "body": {
                "model": "gpt-5",
                "response_format": { 
                    "type": "json_schema",
                    "json_schema": {
                        "name": "paper_extraction",
                        "strict": True,
                        "schema": schema
                    }
                },
                "messages": [
                    {
                        "role": "system",
                        "content": "Extract structured clinical study data from the provided abstracts. Return only valid JSON conforming exactly to the schema, using null for missing values. Convert written numbers into digits where applicable (e.g., forty-six to 46)"
                    },
                    {
                        "role": "user", 
                        "content": json.dumps(papers)
                    }
                ],
                "max_completion_tokens": 20000
            }
        }
        
        f.write(json.dumps(json_object) + '\n')

print(f"Created {(len(df) + batch_size - 1) // batch_size} batched requests")

batch_input_file = client.files.create(
    file=open("batch_requests.jsonl", "rb"),
    purpose="batch"
)

print(f"Uploaded file: {batch_input_file.filename}")
print(f"File ID: {batch_input_file.id}")

batch = client.batches.create(
    input_file_id=batch_input_file.id,
    endpoint="/v1/chat/completions",
    completion_window="24h",
    metadata={"description": "batch"}
)

print("Batch created successfully!")
print(f"Batch ID: {batch.id}")
print(f"Status: {batch.status}")
print(f"Created at: {batch.created_at}")