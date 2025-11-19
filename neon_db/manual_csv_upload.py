#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
from sqlalchemy import create_engine, text
import logging
from datetime import datetime

# Database configuration
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require')

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('manual_upload.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def upload_csv_to_neon(csv_path: str, schema: str, table: str, mode: str = 'append'):
    """Upload CSV file to Neon database"""
    logger = setup_logging()
    
    try:
        # Validate file exists
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        
        # Create database connection
        engine = create_engine(DATABASE_URL)
        logger.info(f"Connected to database")
        
        # Check if table exists
        with engine.connect() as conn:
            table_exists = conn.execute(text("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = :schema 
                    AND table_name = :table
                );
            """), {"schema": schema, "table": table}).scalar()
        
        logger.info(f"Table {schema}.{table} exists: {table_exists}")
        
        # Read and process CSV
        logger.info(f"Reading CSV: {csv_path}")
        
        # For large files, process in chunks
        file_size = os.path.getsize(csv_path)
        chunk_size = 10000 if file_size > 50 * 1024 * 1024 else None  # 50MB threshold
        
        total_rows = 0
        start_time = datetime.now()
        
        if chunk_size:
            logger.info(f"Processing large file in chunks of {chunk_size}")
            for i, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunk_size)):
                # Clean column names
                chunk.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                               for col in chunk.columns]
                
                # Handle null values
                chunk = chunk.where(pd.notnull(chunk), None)
                
                # Upload chunk
                if_exists = mode if i == 0 else 'append'
                chunk.to_sql(table, engine, schema=schema, if_exists=if_exists, 
                           index=False, method='multi')
                
                total_rows += len(chunk)
                logger.info(f"Processed chunk {i+1}, rows so far: {total_rows}")
        else:
            # Process entire file at once
            df = pd.read_csv(csv_path)
            logger.info(f"CSV shape: {df.shape}")
            
            # Clean column names
            df.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                         for col in df.columns]
            
            # Handle null values
            df = df.where(pd.notnull(df), None)
            
            # Upload to database
            df.to_sql(table, engine, schema=schema, if_exists=mode, 
                     index=False, method='multi')
            total_rows = len(df)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        logger.info(f"Upload completed successfully!")
        logger.info(f"Total rows uploaded: {total_rows}")
        logger.info(f"Duration: {duration:.2f} seconds")
        logger.info(f"Target: {schema}.{table}")
        
        return True
        
    except Exception as e:
        logger.error(f"Upload failed: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Upload CSV to Neon Database')
    parser.add_argument('csv_file', help='Path to CSV file')
    parser.add_argument('--schema', default='jim_data', help='Database schema (default: jim_data)')
    parser.add_argument('--table', required=True, help='Table name')
    parser.add_argument('--mode', choices=['append', 'replace', 'fail'], 
                       default='append', help='Upload mode (default: append)')
    
    args = parser.parse_args()
    
    print(f"Uploading {args.csv_file} to {args.schema}.{args.table}")
    print(f"Mode: {args.mode}")
    
    # Confirm before proceeding
    response = input("Continue? (y/N): ")
    if response.lower() != 'y':
        print("Cancelled")
        sys.exit(0)
    
    success = upload_csv_to_neon(args.csv_file, args.schema, args.table, args.mode)
    
    if success:
        print("✅ Upload successful!")
        sys.exit(0)
    else:
        print("❌ Upload failed!")
        sys.exit(1)

if __name__ == '__main__':
    main()