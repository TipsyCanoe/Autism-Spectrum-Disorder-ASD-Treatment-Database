#!/usr/bin/env python3
import os
import pandas as pd
from sqlalchemy import create_engine
import logging
from datetime import datetime

# Database configuration
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require&channel_binding=require')

def interactive_upload():
    """Interactive CSV upload script"""
    print("Neon Database CSV Uploader")
    print("=" * 40)
    
    # Get input from user with PubmedAPIFiles/ prepended
    csv_filename = input("Enter CSV filename (in PubmedAPIFiles/): ").strip().strip('"').strip("'")
    csv_path = f"PubmedAPIFiles/{csv_filename}"
    
    # Validate file
    if not os.path.exists(csv_path):
        print(f"File not found: {csv_path}")
        return False
    
    if not csv_path.lower().endswith('.csv'):
        print("File must be a CSV")
        return False
    
    # Get database details with defaults
    schema = input("Enter schema name [jim_data]: ").strip() or 'jim_data'
    table = input("Enter table name [data_embedded]: ").strip() or 'data_embedded'
    mode = input("Upload mode (append/replace/upsert) [append]: ").strip().lower() or 'append'
    
    if mode not in ['append', 'replace', 'upsert']:
        print("Mode must be 'append', 'replace', or 'upsert'")
        return False
    
    # Preview CSV
    try:
        df_preview = pd.read_csv(csv_path, nrows=5, lineterminator='\n', encoding='utf-8')
        # Strip whitespace from column names
        df_preview.columns = df_preview.columns.str.strip()
        print(f"\nCSV Preview ({df_preview.shape[0]} rows shown):")
        print(df_preview.to_string(index=False))
        print(f"\nColumns: {list(df_preview.columns)}")
        
        # Get file size
        file_size = os.path.getsize(csv_path)
        print(f"File size: {file_size / 1024 / 1024:.2f} MB")
        
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return False
    
    # Confirm upload
    print(f"\nUpload Summary:")
    print(f"   File: {csv_path}")
    print(f"   Target: {schema}.{table}")
    print(f"   Mode: {mode}")
    
    confirm = input("\nProceed with upload? (y/N): ").lower()
    if confirm != 'y':
        print("Upload cancelled")
        return False
    
    # Perform upload
    try:
        print("\nStarting upload...")
        engine = create_engine(DATABASE_URL)
        
        start_time = datetime.now()
        
        # Read full CSV
        df = pd.read_csv(csv_path, lineterminator='\n', encoding='utf-8', keep_default_na=True, na_values=['', 'NaN', 'nan', 'None'])
        # Strip whitespace from column names
        df.columns = df.columns.str.strip()
        print(f"   Loaded {len(df)} rows, {len(df.columns)} columns")
        
        # Clean column names
        original_columns = df.columns.tolist()
        df.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                     for col in df.columns]
        
        # Convert data types to match database schema
        # pmid and sample_size should be integers
        if 'pmid' in df.columns:
            df['pmid'] = pd.to_numeric(df['pmid'], errors='coerce').astype('Int64')
        if 'sample_size' in df.columns:
            df['sample_size'] = pd.to_numeric(df['sample_size'], errors='coerce').astype('Int64')
        
        # Handle null values - replace pandas NaN/None with Python None
        df = df.where(pd.notnull(df), None)
        
        # Special handling for embedding column - ensure NaN becomes None (NULL in database)
        if 'embedding' in df.columns:
            df['embedding'] = df['embedding'].apply(lambda x: None if pd.isna(x) or x == 'NaN' or x == '' else x)
        
        # Upload to database
        if mode == 'upsert':
            # Use custom upsert logic with ON CONFLICT
            from sqlalchemy import text
            
            # Get columns
            columns = df.columns.tolist()
            cols_str = ', '.join(columns)
            placeholders = ', '.join([f':{col}' for col in columns])
            
            # Build update clause (update all columns except pmid)
            update_cols = [col for col in columns if col != 'pmid']
            update_str = ', '.join([f'{col} = EXCLUDED.{col}' for col in update_cols])
            
            # Create upsert SQL
            upsert_sql = text(f"""
                INSERT INTO {schema}.{table} ({cols_str})
                VALUES ({placeholders})
                ON CONFLICT (pmid) DO UPDATE SET {update_str}
            """)
            
            # Execute batch upsert
            with engine.begin() as conn:
                for _, row in df.iterrows():
                    conn.execute(upsert_sql, row.to_dict())
        else:
            # Use standard pandas to_sql for append/replace
            df.to_sql(table, engine, schema=schema, if_exists=mode, 
                     index=False, method='multi')
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"Upload successful!")
        print(f"   Rows uploaded: {len(df)}")
        print(f"   Duration: {duration:.2f} seconds")
        print(f"   Table: {schema}.{table}")
        
        return True
        
    except Exception as e:
        print(f"Upload failed: {e}")
        return False

if __name__ == '__main__':
    success = interactive_upload()
    if not success:
        exit(1)