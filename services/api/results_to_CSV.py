import json
import pandas as pd

all_results = []
errors = []

with open("outputs.jsonl", "r") as f:
    for line_num, line in enumerate(f, 1):
        try:
            result = json.loads(line)
            
            if result['response']['status_code'] == 200:
                pmid = result['custom_id']
                
                content = result['response']['body']['choices'][0]['message']['content']
                paper_data = json.loads(content)
                
                paper_data['pmid'] = pmid
                
                all_results.append(paper_data)
            else:
                errors.append({
                    'custom_id': result['custom_id'],
                    'status_code': result['response']['status_code'],
                    'error': result.get('error', 'Unknown error')
                })
                
        except Exception as e:
            errors.append({
                'line': line_num,
                'error': str(e),
                'content': line[:100] 
            })

df_results = pd.DataFrame(all_results)
df_results.to_csv('extracted_papers.csv', index=False)
print(f"Saved {len(df_results)} papers to extracted_papers.csv")

if errors:
    df_errors = pd.DataFrame(errors)
    df_errors.to_csv('extraction_errors.csv', index=False)
    print(f"{len(errors)} errors saved to extraction_errors.csv")