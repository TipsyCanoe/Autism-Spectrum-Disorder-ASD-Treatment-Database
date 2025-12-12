import pandas as pd
import json
import time
from datetime import datetime, timedelta
from openai import OpenAI
from transformers import AutoTokenizer
from sentence_transformers import SentenceTransformer, models
from tqdm import tqdm
from dotenv import load_dotenv
import os

class LLMPipeline:
    def __init__(self, input_csv='grabbed_papers.csv', batch_size=5):
        self.client = OpenAI()
        self.input_csv = input_csv
        self.batch_size = batch_size
        self.batch_id = None
        self.schema = {
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
                            "pmid", "sample_size", "treatment_name", "study_type",
                            "study_duration", "male_female_ratio", "age_range",
                            "treatment_dose_range", "primary_outcome_area",
                            "primary_outcome_measures", "results_primary_measure",
                            "secondary_outcome_area", "secondary_outcome_measures",
                            "results_secondary_measures", "tolerability_side_effects",
                            "safety", "drop_out_rate", "ethnicity_percentages"
                        ],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["results"],
            "additionalProperties": False
        }
    
    def step1_create_batch(self):
        """Create batch requests from input CSV"""
        print("\n" + "="*50)
        print("STEP 1: Creating batch requests")
        print("="*50)
        
        df = pd.read_csv(self.input_csv)
        print(f"Loaded {len(df)} papers from {self.input_csv}")
        
        output_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/batch_requests.jsonl')
        with open(output_path, 'w', encoding='utf-8') as f:
            for i in range(0, len(df), self.batch_size):
                batch = df.iloc[i:i+self.batch_size]
                
                papers = []
                for _, row in batch.iterrows():
                    papers.append({
                        "pmid": str(row['pmid']) if pd.notna(row['pmid']) else 'N/A',
                        "abstract": row['abstract'] if pd.notna(row['abstract']) else 'N/A'
                    })
                
                json_object = {
                    "custom_id": f"batch_{i//self.batch_size}",
                    "method": "POST",
                    "url": "/v1/chat/completions",
                    "body": {
                        "model": "gpt-5",
                        "response_format": { 
                            "type": "json_schema",
                            "json_schema": {
                                "name": "paper_extraction",
                                "strict": True,
                                "schema": self.schema
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
        
        num_batches = (len(df) + self.batch_size - 1) // self.batch_size
        print(f"Created {num_batches} batched requests in batch_requests.jsonl")
        return True
    
    def step2_send_to_api(self):
        """Upload batch file and create batch job"""
        print("\n" + "="*50)
        print("STEP 2: Sending batch to API")
        print("="*50)
        
        input_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/batch_requests.jsonl')
        batch_input_file = self.client.files.create(
            file=open(input_path, "rb"),
            purpose="batch"
        )
        
        print(f"Uploaded file: {batch_input_file.filename}")
        print(f"File ID: {batch_input_file.id}")
        
        batch = self.client.batches.create(
            input_file_id=batch_input_file.id,
            endpoint="/v1/chat/completions",
            completion_window="24h",
            metadata={"description": f"Research pipeline batch {datetime.now().isoformat()}"}
        )
        
        self.batch_id = batch.id
        
        print(f" Batch created successfully!")
        print(f"  Batch ID: {batch.id}")
        print(f"  Status: {batch.status}")
        print(f"  Created at: {batch.created_at}")
        
        # Save batch ID to file for recovery
        info_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/batch_info.json')
        with open(info_path, 'w') as f:
            json.dump({
                'batch_id': batch.id,
                'created_at': str(datetime.now()),
                'file_id': batch_input_file.id
            }, f, indent=2)
        
        print(f"Batch info saved to batch_info.json")
        return batch.id
    
    def step3_wait_and_check(self, batch_id=None, check_interval=300):
        """Wait for batch to complete and check status periodically"""
        print("\n" + "="*50)
        print("STEP 3: Monitoring batch progress")
        print("="*50)
        
        if batch_id is None:
            batch_id = self.batch_id
        
        if batch_id is None:
            # Try to load from file
            if os.path.exists('batch_info.json'):
                with open('batch_info.json', 'r') as f:
                    info = json.load(f)
                    batch_id = info['batch_id']
                    print(f"Loaded batch ID from file: {batch_id}")
            else:
                raise ValueError("No batch_id provided and no batch_info.json found")
        
        print(f"Monitoring batch: {batch_id}")
        print(f"Check interval: {check_interval} seconds ({check_interval/60:.1f} minutes)")
        
        while True:
            batch = self.client.batches.retrieve(batch_id)
            
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]")
            print(f"  Status: {batch.status}")
            
            if hasattr(batch, 'request_counts'):
                counts = batch.request_counts
                print(f"  Completed: {counts.completed}/{counts.total}")
                print(f"  Failed: {counts.failed}")
            
            if batch.status == "completed":
                print("\nBatch processing completed!")
                return batch
            elif batch.status == "failed" or batch.status == "expired" or batch.status == "cancelled":
                print(f"\nBatch {batch.status}!")
                if hasattr(batch, 'errors'):
                    print(f"  Errors: {batch.errors}")
                raise Exception(f"Batch processing {batch.status}")
            
            print(f"  Waiting {check_interval} seconds before next check...")
            time.sleep(check_interval)
    
    def step4_download_results(self, batch_id=None):
        """Download and save batch results"""
        print("\n" + "="*50)
        print("STEP 4: Downloading results")
        print("="*50)
        
        if batch_id is None:
            batch_id = self.batch_id
        
        if batch_id is None:
            info_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/batch_info.json')
            if os.path.exists(info_path):
                with open(info_path, 'r') as f:
                    info = json.load(f)
                    batch_id = info['batch_id']
        
        batch = self.client.batches.retrieve(batch_id)
        
        if batch.output_file_id is None:
            raise Exception("No output file available")
        
        print(f"Downloading output file: {batch.output_file_id}")
        
        file_response = self.client.files.content(batch.output_file_id)
        
        output_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/outputs.jsonl')
        with open(output_path, 'wb') as f:
            f.write(file_response.content)
        
        print(f"Results saved to {output_path}")
        return output_path
    
    def step5_convert_to_csv(self, input_file=None):
        """Convert JSONL results to CSV"""
        print("\n" + "="*50)
        print("STEP 5: Converting results to CSV")
        print("="*50)
        
        if input_file is None:
            input_file = os.path.join(os.path.dirname(__file__), '../../data/embeddings/outputs.jsonl')
        
        all_results = []
        errors = []
        
        with open(input_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                try:
                    result = json.loads(line)
                    
                    if result['response']['status_code'] == 200:
                        content = result['response']['body']['choices'][0]['message']['content']
                        data = json.loads(content)
                        
                        # Handle both single results and batched results
                        if 'results' in data:
                            # Batched format
                            for item in data['results']:
                                item['batch_id'] = result['custom_id']
                                all_results.append(item)
                        else:
                            # Single format
                            data['batch_id'] = result['custom_id']
                            all_results.append(data)
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
        output_csv = os.path.join(os.path.dirname(__file__), '../../data/embeddings/extracted_papers.csv')
        df_results.to_csv(output_csv, index=False)
        print(f"Saved {len(df_results)} papers to {output_csv}")
        
        if errors:
            df_errors = pd.DataFrame(errors)
            error_csv = os.path.join(os.path.dirname(__file__), '../../data/embeddings/extraction_errors.csv')
            df_errors.to_csv(error_csv, index=False)
            print(f"{len(errors)} errors saved to {error_csv}")
        
        return output_csv
    
    def step6_create_vectors(self, input_csv=None, 
                            output_csv=None,
                            model_name='Charangan/MedBERT'):
        """Generate embeddings for abstracts"""
        print("\n" + "="*50)
        print("STEP 6: Creating vector embeddings")
        print("="*50)
        
        if input_csv is None:
            input_csv = os.path.join(os.path.dirname(__file__), '../../data/embeddings/extracted_papers.csv')
        if output_csv is None:
            output_csv = os.path.join(os.path.dirname(__file__), '../../data/embeddings/LLM_with_embeddings.csv')
        
        df = pd.read_csv(input_csv, encoding="utf-8", engine="python")
        print(f"Loaded {len(df)} papers from {input_csv}")
        
        if "abstract" not in df.columns:
            raise ValueError("CSV must contain an 'abstract' column.")
        
        print(f"Loading model: {model_name}")
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        word_embedding_model = models.Transformer(model_name)
        pooling_model = models.Pooling(
            word_embedding_model.get_word_embedding_dimension(),
            pooling_mode_mean_tokens=True,
            pooling_mode_cls_token=False,
            pooling_mode_max_tokens=False
        )
        sentence_model = SentenceTransformer(modules=[word_embedding_model, pooling_model])
        
        print("Encoding abstracts...")
        embeddings = []
        for text in tqdm(df["abstract"].fillna(""), desc="Encoding abstracts"):
            if isinstance(text, str) and text.strip():
                emb = sentence_model.encode(text)
                embeddings.append(emb.tolist())
            else:
                embeddings.append(None)
        
        df["vector"] = embeddings
        df["AI"] = True
        df["vector"] = df["vector"].apply(lambda v: str(v) if v is not None else "")
        
        df.to_csv(output_csv, index=False, encoding="utf-8")
        print(f"Saved embeddings to {output_csv}")
        
        return output_csv
    
    def run_full_pipeline(self, wait_for_completion=True, check_interval=300, wait_24h=False):
        """Execute the complete pipeline
        
        Args:
            wait_for_completion: Whether to wait for batch completion
            check_interval: How often to check status (seconds) when actively monitoring
            wait_24h: If True, sleep for 24 hours before checking (ignores check_interval)
        """
        print("\n" + "="*70)
        print("  RESEARCH PAPER PROCESSING PIPELINE")
        print("="*70)
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        start_time = time.time()
        
        try:
            # Step 1: Create batch requests
            self.step1_create_batch()
            
            # Step 2: Send to API
            batch_id = self.step2_send_to_api()
            
            # Step 3: Wait for completion (or skip if not waiting)
            if wait_for_completion:
                wait_until = datetime.now() + timedelta(hours=24)
                print("\n" + "="*50)
                print("SLEEPING FOR 24 HOURS")
                print("="*50)
                print(f"Current time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"Will check at: {wait_until.strftime('%Y-%m-%d %H:%M:%S')}")
                print("\nYou can safely close this terminal and resume later with:")
                print(f"  pipeline.resume_from_batch('{batch_id}')")
                
                time.sleep(24 * 60 * 60)  # Sleep for 24 hours
                
                print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Woke up after 24 hours")
                print("Checking batch status...")
                
                batch = self.client.batches.retrieve(batch_id)
                if batch.status != "completed":
                    print(f"Status: {batch.status} - Continuing to monitor...")
                    self.step3_wait_and_check(batch_id, check_interval)
                
                # Step 4: Download results
                self.step4_download_results(batch_id)
                
                # Step 5: Convert to CSV
                self.step5_convert_to_csv()
                
                # Step 6: Create vectors
                self.step6_create_vectors()
                
                elapsed = time.time() - start_time
                print("\n" + "="*70)
                print("  PIPELINE COMPLETED SUCCESSFULLY!")
                print("="*70)
                print(f"Total time: {elapsed/3600:.2f} hours ({elapsed/60:.1f} minutes)")
                print(f"Final output: LLM_with_embeddings.csv")
            else:
                print("\n" + "="*70)
                print("  BATCH SUBMITTED - Waiting skipped")
                print("="*70)
                print(f"Batch ID: {batch_id}")
                print("\nTo resume later, run:")
                print(f"  pipeline.resume_from_batch('{batch_id}')")
                
        except Exception as e:
            print("\n" + "="*70)
            print("  PIPELINE FAILED")
            print("="*70)
            print(f"Error: {str(e)}")
            raise
    
    def resume_from_batch(self, batch_id=None, check_interval=300):
        """Resume pipeline from an existing batch"""
        print("\n" + "="*70)
        print("  RESUMING PIPELINE FROM EXISTING BATCH")
        print("="*70)
        
        if batch_id is None:
            info_path = os.path.join(os.path.dirname(__file__), '../../data/embeddings/batch_info.json')
            if os.path.exists(info_path):
                with open(info_path, 'r') as f:
                    info = json.load(f)
                    batch_id = info['batch_id']
                    print(f"Loaded batch ID from file: {batch_id}")
            else:
                raise ValueError("No batch_id provided and no batch_info.json found")
        
        # Check current status
        batch = self.client.batches.retrieve(batch_id)
        print(f"Current status: {batch.status}")
        
        if batch.status == "completed":
            print("Batch already completed, proceeding to download...")
            self.step4_download_results(batch_id)
            self.step5_convert_to_csv()
            self.step6_create_vectors()
            print("\nPipeline resumed and completed successfully!")
        else:
            print("Waiting for batch to complete...")
            self.step3_wait_and_check(batch_id, check_interval)
            self.step4_download_results(batch_id)
            self.step5_convert_to_csv()
            self.step6_create_vectors()
            print("\nPipeline resumed and completed successfully!")


if __name__ == "__main__":
    load_dotenv(override=True)
    last_pull_date = os.getenv('LAST_PULL_DATE')
    
    # Define paths
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_RAW_DIR = os.path.join(SCRIPT_DIR, "../../data/raw")
    
    input_file = os.path.join(DATA_RAW_DIR, f"pubmed_combined_{last_pull_date}.csv")
    
    pipeline = LLMPipeline(
        input_csv=input_file,
        batch_size=5
    )
    
    pipeline.run_full_pipeline(wait_for_completion=True)