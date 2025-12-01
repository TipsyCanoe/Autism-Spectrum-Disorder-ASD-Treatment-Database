#!/usr/bin/env python3
"""
Automated upload script for LLM_with_embeddings.csv
Hardcoded to upload to jim_data.data_embedded in append mode
"""
import os
import pandas as pd
from sqlalchemy import create_engine
from datetime import datetime

# Database configuration
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require&channel_binding=require')

# Hardcoded configuration
CSV_FILENAME = 'LLM_with_embeddings.csv'
CSV_PATH = f'PubmedAPIFiles/{CSV_FILENAME}'
SCHEMA = 'jim_data'
TABLE = 'data_embedded'
MODE = 'append'

def upload_llm_embeddings():
    """Upload LLM_with_embeddings.csv to jim_data.data_embedded"""
    print("=" * 50)
    print("LLM Embeddings Auto-Upload to Neon Database")
    print("=" * 50)
    print(f"\nConfiguration:")
    print(f"  File: {CSV_PATH}")
    print(f"  Target: {SCHEMA}.{TABLE}")
    print(f"  Mode: {MODE}")
    print("-" * 50)
    
    # Validate file exists
    if not os.path.exists(CSV_PATH):
        print(f"Error: File not found: {CSV_PATH}")
        print(f"   Please ensure {CSV_FILENAME} exists in the PubmedAPIFiles/ directory")
        return False
    
    if not CSV_PATH.lower().endswith('.csv'):
        print(f"Error: File must be a CSV")
        return False
    
    # Preview CSV
    try:
        df_preview = pd.read_csv(CSV_PATH, nrows=5)
        print(f"\nCSV Preview ({df_preview.shape[0]} rows shown):")
        print(df_preview.to_string(index=False, max_colwidth=50))
        print(f"\nColumns: {list(df_preview.columns)}")
        
        # Get file size
        file_size = os.path.getsize(CSV_PATH)
        print(f"File size: {file_size / 1024 / 1024:.2f} MB")
        
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return False
    
    # Perform upload
    try:
        print(f"\nStarting upload...")
        engine = create_engine(DATABASE_URL)
        
        start_time = datetime.now()
        
        # Read full CSV
        df = pd.read_csv(CSV_PATH)
        print(f"   Loaded {len(df)} rows, {len(df.columns)} columns")
        
        # Clean column names
        df.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                     for col in df.columns]
        
        # Handle null values
        df = df.where(pd.notnull(df), None)
        
        # Upload to database
        print(f"   Uploading to {SCHEMA}.{TABLE}...")
        df.to_sql(TABLE, engine, schema=SCHEMA, if_exists=MODE, 
                 index=False, method='multi')
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"\nUpload successful!")
        print(f"   Rows uploaded: {len(df)}")
        print(f"   Duration: {duration:.2f} seconds")
        print(f"   Table: {SCHEMA}.{TABLE}")
        print(f"   Timestamp: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 50)
        
        return True
        
    except Exception as e:
        print(f"\nUpload failed: {e}")
        print("=" * 50)
        return False

if __name__ == '__main__':
    success = upload_llm_embeddings()
    if not success:
        exit(1)
