from flask import Flask, request, jsonify
from flask_cors import CORS
import psycopg2
from psycopg2.extras import execute_values
import numpy as np
from transformers import AutoModel, AutoTokenizer
from sentence_transformers import SentenceTransformer, models
import os
from dotenv import load_dotenv
import ssl
import re


app = Flask(__name__)
CORS(app) 

NEON_DATABASE_URL = os.getenv("DATABASE_URL")

from transformers import AutoModel, AutoTokenizer
from sentence_transformers import SentenceTransformer, models

medbert_model = AutoModel.from_pretrained("Charangan/MedBERT")
tokenizer = AutoTokenizer.from_pretrained("Charangan/MedBERT")

word_embedding_model = models.Transformer("Charangan/MedBERT")
pooling_model = models.Pooling(
    word_embedding_model.get_word_embedding_dimension(),
    pooling_mode_mean_tokens=True,
    pooling_mode_cls_token=False,
    pooling_mode_max_tokens=False
)

sentence_model = SentenceTransformer(modules=[word_embedding_model, pooling_model])

available_filters = {
    "age": ["0-5", "6-12", "13-17", "18-25", "26-64", "65+"],
    "symptom": [
        "irritability", "adhd", "hyperactivity", "social", 
        "attention-hyperactivity",
        "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance", 
        "anxiety-reactivity"
    ],
    "gender": ["male", "female", "nonbinary"],
    "medication": ["aripiprazole", "citalopram"]
}

def extract_age_from_range(age_range_str):
    """Extract age information from age range string like '13-17 / 15' or '10-17 / 14.2'"""
    if not age_range_str:
        return None
    
    # Extract the range part before the '/'
    range_part = age_range_str.split('/')[0].strip()
    
    # Extract numbers from the range
    numbers = re.findall(r'\d+', range_part)
    if len(numbers) >= 2:
        min_age, max_age = int(numbers[0]), int(numbers[1])
        return min_age, max_age
    return None

def matches_age_filter(study_age_range, filter_age_ranges):
    """Check if study age range overlaps with any of the filter age ranges"""
    study_ages = extract_age_from_range(study_age_range)
    if not study_ages:
        return False
    
    study_min, study_max = study_ages
    
    for filter_range in filter_age_ranges:
        if filter_range == "0-5":
            filter_min, filter_max = 0, 5
        elif filter_range == "6-12":
            filter_min, filter_max = 6, 12
        elif filter_range == "13-17":
            filter_min, filter_max = 13, 17
        elif filter_range == "18-25":
            filter_min, filter_max = 18, 25
        elif filter_range == "26-64":
            filter_min, filter_max = 26, 64
        elif filter_range == "65+":
            filter_min, filter_max = 65, 150
        else:
            continue
        
        # Check for overlap
        if study_min <= filter_max and study_max >= filter_min:
            return True
    
    return False

def matches_gender_filter(study_mf_ratio, filter_genders):
    """Check if study includes the filtered genders based on M:F ratio"""
    if not study_mf_ratio:
        return False
    
    if "male" in filter_genders or "female" in filter_genders:
        return True
    
    return False

def matches_symptom_filter(study, filter_symptoms):
    """Check if study relates to any of the filtered symptoms"""
    searchable_text = " ".join([
        study.get("Primary Outcome Area", ""),
        study.get("Secondary Outcome Area", ""),
        study.get("Study Title", ""),
        study.get("Results: Primary measure", ""),
        study.get("Results: Secondary Measures", "")
    ]).lower()
    
    for symptom in filter_symptoms:
        symptom_terms = {
            "irritability": ["irritability", "irritable"],
            "adhd": ["adhd", "attention deficit", "hyperactivity disorder"],
            "hyperactivity": ["hyperactivity", "hyperactive"],
            "social": ["social", "autism", "asd"],
            "attention-hyperactivity": ["attention", "hyperactivity"],
            "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance": ["lethargy", "withdrawal", "stereotypy"],
            "anxiety-reactivity": ["anxiety", "anxious", "reactivity"]
        }
        
        terms = symptom_terms.get(symptom, [symptom])
        if any(term in searchable_text for term in terms):
            return True
    
    return False

def get_db_connection():
    """Create a connection to the Neon PostgreSQL database"""
    # SSL context for secure connection
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE
    
    # Option 2: Using connection string
    if NEON_DATABASE_URL:
        conn = psycopg2.connect(NEON_DATABASE_URL, sslmode='require', sslcontext=ssl_context)
    else:
        raise ValueError("Database connection information not provided")
        
    return conn

def matches_text_query(study, query):
    """Check if study matches the text query"""
    if not query:
        return True
    
    # Combine all searchable fields
    searchable_text = " ".join([
        study.get("Study Title", ""),
        study.get("Primary Outcome Area", ""),
        study.get("Secondary Outcome Area", ""),
        study.get("Results: Primary measure", ""),
        study.get("Results: Secondary Measures", ""),
        study.get("Tolerability/Side Effects", ""),
        study.get("Safety", "")
    ]).lower()
    
    query_terms = query.lower().split()
    return all(term in searchable_text for term in query_terms)

@app.route('/api/filters', methods=['GET'])
def get_filters():
    """Return available filter options"""
    return jsonify(available_filters)

@app.route('/api/search', methods=['GET'])
def search():
    """Search resources based on query and filters using vector similarity"""
    query = request.args.get('query', '').lower()
    selected_filters = request.args.getlist('filters')
    limit = int(request.args.get('limit', 10))
    
    # Parse filters
    parsed_filters = {}
    for filter_str in selected_filters:
        if ':' in filter_str:
            category, value = filter_str.split(':', 1)
            if category not in parsed_filters:
                parsed_filters[category] = []
            parsed_filters[category].append(value)
    
    conn = get_db_connection()
    cursor = conn.cursor()
    
    results = []
    
    if query:
        # Generate embedding for the query
        query_embedding = sentence_model.encode(query)
        
        # Build the SQL query with filters
        sql = """
        SELECT id, title, type, publish_date, author, description, tags, url, age, symptom, gender,
               embedding <=> %s AS distance
        FROM resources
        WHERE 1=1
        """
        params = [query_embedding.tolist()]
        
        # Add filter conditions
        for category, values in parsed_filters.items():
            if values:
                placeholders = ', '.join(['%s'] * len(values))
                sql += f" AND {category} IN ({placeholders})"
                params.extend(values)
        
        # Order by vector similarity (cosine distance)
        sql += " ORDER BY distance ASC LIMIT %s"
        params.append(limit)
        
        cursor.execute(sql, params)
        
        # Process results
        for row in cursor.fetchall():
            results.append({
                "id": row[0],
                "title": row[1],
                "type": row[2],
                "publishDate": row[3].strftime('%Y-%m-%d'),
                "author": row[4],
                "description": row[5],
                "tags": row[6],
                "url": row[7],
                "age": row[8],
                "symptom": row[9],
                "gender": row[10],
                "similarity_score": 1 - float(row[11])  # Convert distance to similarity score
            })
    else:
        # If no query, just apply filters
        sql = """
        SELECT id, title, type, publish_date, author, description, tags, url, age, symptom, gender
        FROM resources
        WHERE 1=1
        """
        params = []
        
        # Add filter conditions
        for category, values in parsed_filters.items():
            if values:
                placeholders = ', '.join(['%s'] * len(values))
                sql += f" AND {category} IN ({placeholders})"
                params.extend(values)
        
        sql += " LIMIT %s"
        params.append(limit)
        
        cursor.execute(sql, params)
        
        # Process results
        for row in cursor.fetchall():
            results.append({
                "id": row[0],
                "title": row[1],
                "type": row[2],
                "publishDate": row[3].strftime('%Y-%m-%d'),
                "author": row[4],
                "description": row[5],
                "tags": row[6],
                "url": row[7],
                "age": row[8],
                "symptom": row[9],
                "gender": row[10]
            })
    
    cursor.close()
    conn.close()
    print(jsonify({"results": results, "total": len(results)}))
    return jsonify({"results": results, "total": len(results)})

if __name__ == '__main__':
    app.run(debug=True, port=5000)