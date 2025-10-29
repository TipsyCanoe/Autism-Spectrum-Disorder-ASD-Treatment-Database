import pandas as pd
import json

df = pd.read_csv('grabbed_papers.csv')

with open('batch_requests.jsonl', 'w', encoding='utf-8') as f:
    for _, row in df.iterrows():
        pmid = row['pmid'] if pd.notna(row['pmid']) else 'N/A'
        title = row['title'] if pd.notna(row['title']) else 'N/A'
        abstract = row['abstract'] if pd.notna(row['abstract']) else 'N/A'
        keywords = row['keywords'] if pd.notna(row['keywords']) else 'N/A'
        
        formatted_string = f"pmid: {pmid}, title: {title}, abstract: {abstract}, keywords: {keywords}"
        
        json_object = {
            "custom_id": str(pmid),
            "method": "POST",
            "url": "/v1/chat/completions",
            "body": {
                "model": "gpt-5",
                "messages": [
                    {
                        "role": "system",
                        "content": "give me a json file with the treatment name, study type, study duration, Male:female ratio, age range, treatment dose range, primary outcome area, primary outcome measures, results: primary measure, secondary outcome area, secondary outcome measures, results:secondary measures, tolerability/side effects, safety, drop out rate, and ethnicity percentages. Please add null if you don't know the answer."
                    },
                    {
                        "role": "user", 
                        "content": formatted_string
                    }
                ],
                "max_completion_tokens": 1000
            }
        }
        
        f.write(json.dumps(json_object) + '\n')

print(f"Created JSONL file with {len(df)} requests")