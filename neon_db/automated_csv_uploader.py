import os
import pandas as pd
import psycopg2
from psycopg2 import sql
from sqlalchemy import create_engine, text
from flask import Flask, request, jsonify
from werkzeug.utils import secure_filename
import logging
from datetime import datetime
import json
from typing import Dict, List, Optional
import tempfile

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('csv_upload.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100MB max file size

# Database configuration - use your existing connection string
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require')

class CSVUploader:
    def __init__(self, database_url: str):
        self.database_url = database_url
        self.engine = create_engine(database_url)
    
    def validate_csv(self, file_path: str) -> Dict:
        """Validate CSV file and return metadata"""
        try:
            # Read first few rows to validate structure
            df = pd.read_csv(file_path, nrows=5)
            
            return {
                'valid': True,
                'columns': df.columns.tolist(),
                'shape': df.shape,
                'dtypes': df.dtypes.to_dict()
            }
        except Exception as e:
            logger.error(f"CSV validation failed: {str(e)}")
            return {'valid': False, 'error': str(e)}
    
    def check_table_exists(self, schema: str, table: str) -> bool:
        """Check if table exists in database"""
        try:
            with self.engine.connect() as conn:
                result = conn.execute(text("""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables 
                        WHERE table_schema = :schema 
                        AND table_name = :table
                    );
                """), {"schema": schema, "table": table})
                return result.scalar()
        except Exception as e:
            logger.error(f"Error checking table existence: {str(e)}")
            return False
    
    def get_table_columns(self, schema: str, table: str) -> List[str]:
        """Get existing table columns"""
        try:
            with self.engine.connect() as conn:
                result = conn.execute(text("""
                    SELECT column_name 
                    FROM information_schema.columns 
                    WHERE table_schema = :schema 
                    AND table_name = :table
                    ORDER BY ordinal_position;
                """), {"schema": schema, "table": table})
                return [row[0] for row in result.fetchall()]
        except Exception as e:
            logger.error(f"Error getting table columns: {str(e)}")
            return []
    
    def create_table_from_csv(self, file_path: str, schema: str, table: str) -> bool:
        """Create table based on CSV structure"""
        try:
            df = pd.read_csv(file_path, nrows=1000)  # Sample for type inference
            
            # Clean column names
            df.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                         for col in df.columns]
            
            # Use pandas to_sql to create table
            df.to_sql(f'{table}_temp', self.engine, schema=schema, 
                     if_exists='replace', index=False)
            
            # Rename temp table to actual table
            with self.engine.connect() as conn:
                conn.execute(text(f'ALTER TABLE {schema}.{table}_temp RENAME TO {table}'))
                conn.commit()
            
            logger.info(f"Created table {schema}.{table}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating table: {str(e)}")
            return False
    
    def upload_csv_data(self, file_path: str, schema: str, table: str, 
                       mode: str = 'append') -> Dict:
        """Upload CSV data to database"""
        try:
            start_time = datetime.now()
            
            # Read CSV in chunks for memory efficiency
            chunk_size = 10000
            chunks_processed = 0
            total_rows = 0
            
            for chunk in pd.read_csv(file_path, chunksize=chunk_size):
                # Clean column names to match database
                chunk.columns = [col.strip().lower().replace(' ', '_').replace('-', '_') 
                               for col in chunk.columns]
                
                # Handle data types and null values
                chunk = chunk.where(pd.notnull(chunk), None)
                
                # Upload chunk
                chunk.to_sql(table, self.engine, schema=schema, 
                           if_exists=mode if chunks_processed == 0 else 'append', 
                           index=False, method='multi')
                
                chunks_processed += 1
                total_rows += len(chunk)
                logger.info(f"Processed chunk {chunks_processed}, total rows: {total_rows}")
            
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            return {
                'success': True,
                'rows_uploaded': total_rows,
                'chunks_processed': chunks_processed,
                'duration_seconds': duration,
                'timestamp': end_time.isoformat()
            }
            
        except Exception as e:
            logger.error(f"Error uploading CSV data: {str(e)}")
            return {'success': False, 'error': str(e)}

# Initialize uploader
uploader = CSVUploader(DATABASE_URL)

@app.route('/upload-csv', methods=['POST'])
def upload_csv():
    """API endpoint to upload CSV files"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        schema = request.form.get('schema', 'jim_data')  # Default schema
        table = request.form.get('table')
        mode = request.form.get('mode', 'append')  # append, replace, or create
        
        if not table:
            return jsonify({'error': 'Table name required'}), 400
        
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if not file.filename.endswith('.csv'):
            return jsonify({'error': 'Only CSV files allowed'}), 400
        
        # Save file temporarily
        filename = secure_filename(file.filename)
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.csv', delete=False) as tmp_file:
            file.save(tmp_file.name)
            temp_path = tmp_file.name
        
        try:
            # Validate CSV
            validation = uploader.validate_csv(temp_path)
            if not validation['valid']:
                return jsonify({'error': f'Invalid CSV: {validation["error"]}'}), 400
            
            # Check if table exists
            table_exists = uploader.check_table_exists(schema, table)
            
            if mode == 'create' or not table_exists:
                if table_exists and mode != 'create':
                    return jsonify({'error': f'Table {schema}.{table} already exists. Use mode="append" or mode="replace"'}), 400
                
                # Create new table
                success = uploader.create_table_from_csv(temp_path, schema, table)
                if not success:
                    return jsonify({'error': 'Failed to create table'}), 500
                
                # Upload data
                result = uploader.upload_csv_data(temp_path, schema, table, 'append')
            else:
                # Table exists, upload data
                upload_mode = 'replace' if mode == 'replace' else 'append'
                result = uploader.upload_csv_data(temp_path, schema, table, upload_mode)
            
            if result['success']:
                response_data = {
                    'message': 'CSV uploaded successfully',
                    'schema': schema,
                    'table': table,
                    'mode': mode,
                    'validation': validation,
                    'upload_result': result
                }
                return jsonify(response_data), 200
            else:
                return jsonify({'error': f'Upload failed: {result["error"]}'}), 500
        
        finally:
            # Clean up temporary file
            os.unlink(temp_path)
    
    except Exception as e:
        logger.error(f"API error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/table-info/<schema>/<table>', methods=['GET'])
def get_table_info(schema, table):
    """Get information about existing table"""
    try:
        if not uploader.check_table_exists(schema, table):
            return jsonify({'error': 'Table does not exist'}), 404
        
        columns = uploader.get_table_columns(schema, table)
        
        # Get row count
        with uploader.engine.connect() as conn:
            result = conn.execute(text(f'SELECT COUNT(*) FROM {schema}.{table}'))
            row_count = result.scalar()
        
        return jsonify({
            'schema': schema,
            'table': table,
            'columns': columns,
            'row_count': row_count
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)