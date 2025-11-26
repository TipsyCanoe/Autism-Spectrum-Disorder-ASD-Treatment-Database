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
    print("🚀 Neon Database CSV Uploader")
    print("=" * 40)
    
    # Get input from user with PubmedAPIFiles/ prepended
    csv_filename = input("Enter CSV filename (in PubmedAPIFiles/): ").strip().strip('"').strip("'")
    csv_path = f"PubmedAPIFiles/{csv_filename}"
    
    # Validate file
    if not os.path.exists(csv_path):
        print(f"❌ File not found: {csv_path}")
        return False
    
    if not csv_path.lower().endswith('.csv'):
        print("❌ File must be a CSV")
        return False
    
    # Get database details with defaults
    schema = input("Enter schema name [jim_data]: ").strip() or 'jim_data'
    table = input("Enter table name [data_embedded]: ").strip() or 'data_embedded'
    mode = input("Upload mode (append/replace) [append]: ").strip().lower() or 'append'
    
    if mode not in ['append', 'replace']:
        print("❌ Mode must be 'append' or 'replace'")
        return False
    
    # Preview CSV
    try:
        df_preview = pd.read_csv(csv_path, nrows=5)
        print(f"\n📊 CSV Preview ({df_preview.shape[0]} rows shown):")
        print(df_preview.to_string(index=False))
        print(f"\nColumns: {list(df_preview.columns)}")
        
        # Get file size
        file_size = os.path.getsize(csv_path)
        print(f"File size: {file_size / 1024 / 1024:.2f} MB")
        
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")
        return False
    
    # Confirm upload
    print(f"\n🎯 Upload Summary:")
    print(f"   File: {csv_path}")
    print(f"   Target: {schema}.{table}")
    print(f"   Mode: {mode}")
    
    confirm = input("\nProceed with upload? (y/N): ").lower()
    if confirm != 'y':
        print("❌ Upload cancelled")
        return False
    
    # Perform upload
    try:
        print("\n🔄 Starting upload...")
        engine = create_engine(DATABASE_URL)
        
        start_time = datetime.now()
        
        # Read full CSV
        df = pd.read_csv(csv_path)
        print(f"   Loaded {len(df)} rows, {len(df.columns)} columns")
        
        # Clean column names
        original_columns = df.columns.tolist()
        df.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                     for col in df.columns]
        
        # Handle null values
        df = df.where(pd.notnull(df), None)
        
        # Upload to database
        df.to_sql(table, engine, schema=schema, if_exists=mode, 
                 index=False, method='multi')
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"✅ Upload successful!")
        print(f"   Rows uploaded: {len(df)}")
        print(f"   Duration: {duration:.2f} seconds")
        print(f"   Table: {schema}.{table}")
        
        return True
        
    except Exception as e:
        print(f"❌ Upload failed: {e}")
        return False

if __name__ == '__main__':
    success = interactive_upload()
    if not success:
        exit(1)