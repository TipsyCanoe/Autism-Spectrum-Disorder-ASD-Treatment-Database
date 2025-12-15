import os
import psycopg2
from dotenv import load_dotenv
from pathlib import Path

# Load env from config directory relative to this script
script_dir = Path(__file__).parent
config_path = script_dir.parent / 'config' / 'production.env'
load_dotenv(config_path)

def get_db_connection():
    connection_url = os.getenv('DATABASE_URL')
    if 'sslmode=' not in connection_url:
        separator = '&' if '?' in connection_url else '?'
        connection_url = f"{connection_url}{separator}sslmode=require"
    conn = psycopg2.connect(connection_url)
    return conn

conn = get_db_connection()
cursor = conn.cursor()

cursor.execute("""
    SELECT column_name, data_type 
    FROM information_schema.columns 
    WHERE table_schema = 'jim_data' 
    AND table_name = 'data_embedded'
    AND column_name = 'ai';
""")

result = cursor.fetchone()
print(f"Column 'ai' type: {result}")

cursor.execute("SELECT DISTINCT ai FROM jim_data.data_embedded;")
rows = cursor.fetchall()
print("Distinct values:", rows)

conn.close()
