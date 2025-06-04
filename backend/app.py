from flask import Flask, request, jsonify
from flask_cors import CORS
import psycopg2
from psycopg2.extras import execute_values
import numpy as np
from transformers import AutoModel, AutoTokenizer
from sentence_transformers import SentenceTransformer, models
from collections import defaultdict
import os
from dotenv import load_dotenv
import ssl
import datetime
import re
import json


app = Flask(__name__)
CORS(app) 

NEON_DATABASE_URL = "postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require"

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
    "medication": []
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
    print("connecting")
    connection_url = "postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require"
    if 'sslmode=' not in connection_url:
        separator = '&' if '?' in connection_url else '?'
        connection_url = f"{connection_url}{separator}sslmode=require"
    
    conn = psycopg2.connect(connection_url)
        
    return conn

def matches_text_query(study, query):
    """Check if study matches the text query"""
    if not query:
        return True
    
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
    print("trying search")
    """Search resources based on query and filters using vector similarity"""
    query = request.args.get('query', '').lower()
    selected_filters = request.args.getlist('filters')
    filters_str = ', '.join(selected_filters)
    limit = int(request.args.get('limit', 30))
    
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
    query_embedding = sentence_model.encode(query + filters_str)

    similarity_threshold = .999
    
    sql = """
        SELECT 
            spv.title,
            spv.pub_date,
            spv.pmid,
            spv.authors,
            spv.url,
            t.treatment_name,
            t.duration,
            t.primary_outcome_area,
            t.primary_outcome_measures,
            spv.embedding <=> %s::vector AS distance,
            spv.abstract
        FROM updated_treatment_data.semantic_paper_search_view spv
        JOIN updated_treatment_data.treatment_pubmed_link tpl ON spv.pmid = tpl.pmid
        JOIN updated_treatment_data.treatments t ON tpl.treatment_id = t.id
        WHERE spv.embedding <=> %s::vector < %s
        ORDER BY distance ASC
        LIMIT %s
    """
    
    params = [query_embedding.tolist(), query_embedding.tolist(), similarity_threshold, limit]
    
    try:
        cursor.execute(sql, params)
        results = cursor.fetchall()
        
        grouped_results = defaultdict(list)
        
        for row in results:
            paper_data = {
                "Abstract": row[10] if row[10] else "N/A",
                "Primary Outcome Area": row[7] if row[7] else "N/A",
                "Primary Outcome Measure": row[8] if row[8] else "N/A",
                "Treatment Duration": row[6] if row[6] else "N/A",
                "Publication Date": row[1].isoformat() if row[1] else "N/A",
                "Author": row[3] if row[3] else "N/A",
                "Study Title": row[0] if row[0] else "N/A",
                "PMID": str(row[2]) if row[2] else "N/A",
                "Full Text URL": row[4] if row[4] else None,
                "Similarity Score": round(1 - row[9], 3) if row[9] is not None else 0,
                "Distance": round(row[9], 3) if row[9] is not None else 1
            }
            
            treatment_name = row[5].lower() if row[5] else "unknown"
            grouped_results[treatment_name].append(paper_data)
        
        final_results = dict(grouped_results)
        final_results = dict(sorted(final_results.items(), 
                                  key=lambda x: len(x[1]), 
                                  reverse=True))
        
        json_output = []
        for treatment_name, studies in final_results.items():
            treatment_entry = {
                "treatment": treatment_name,
                "studies": studies
            }
            json_output.append(treatment_entry)
        print(json.dumps(json_output, indent=4, sort_keys=False))
        return json.dumps(json_output, indent=4, sort_keys=False)
        
    except Exception as e:
        print(f"Error executing vector search: {e}")
        return {}
    
if __name__ == '__main__':
    app.run(debug=True, port=5000)