from flask import Flask, request, jsonify
from flask_cors import CORS
import psycopg2
from psycopg2.extras import execute_values
import numpy as np
from collections import defaultdict
import os
from dotenv import load_dotenv
import ssl
import datetime
import re
import json
import functools  # For the caching decorator

# Simple cache implementation
class SimpleCache:
    def __init__(self, max_size=100):
        self.cache = {}
        self.max_size = max(1, max_size)  # Ensure max_size is at least 1
    
    def get(self, key):
        return self.cache.get(key)
    
    def set(self, key, value):
        # Basic LRU: if cache is full, remove first item (simplistic approach)
        if len(self.cache) >= self.max_size:
            # Remove oldest item (first key)
            if self.cache:  # Check if cache is not empty
                oldest_key = next(iter(self.cache))
                del self.cache[oldest_key]
        self.cache[key] = value

# Load environment variables
load_dotenv()

# Allow disabling heavy model loading for CI/tests via env var
DISABLE_MODEL_LOADING = os.getenv("DISABLE_MODEL_LOADING", "").lower() in ("1", "true", "yes")

# Create cache instances
search_cache = SimpleCache(max_size=75)  # Cache for search results
filter_cache = SimpleCache(max_size=10)   # Cache for filter options

app = Flask(__name__)
CORS(app) 

# Get database URL from environment variable with fallback
NEON_DATABASE_URL = os.getenv('DATABASE_URL', "postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require")

if DISABLE_MODEL_LOADING:
    class _StubSentenceModel:
        def encode(self, text):
            # Return a small fixed-size vector to satisfy downstream code
            import numpy as _np
            return _np.zeros(4, dtype=float)

    sentence_model = _StubSentenceModel()
else:
    # Import heavy ML libraries only when model loading is enabled
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
    # Try to get from cache first
    cached_filters = filter_cache.get('all_filters')
    if cached_filters:
        print("Returning filters from cache")
        return cached_filters
    
    # Add code to load medication options if they're empty
    if not available_filters["medication"]:
        try:
            conn = get_db_connection()
            cursor = conn.cursor()
            # Get distinct treatment names
            cursor.execute("""
                SELECT DISTINCT treatment_name 
                FROM updated_treatment_data.treatments 
                WHERE treatment_name IS NOT NULL AND treatment_name <> '' 
                ORDER BY treatment_name ASC
            """)
            medications = [row[0].lower() for row in cursor.fetchall()]
            available_filters["medication"] = medications
            print(f"Loaded {len(medications)} medications for filter")
        except Exception as e:
            print(f"Error loading medications for filter: {e}")
        finally:
            if cursor:
                cursor.close()
            if conn:
                conn.close()
    
    # Cache the result
    result = jsonify(available_filters)
    filter_cache.set('all_filters', result)
    return result

@app.route('/api/search', methods=['GET'])
def search():
    print("trying search")
    query = request.args.get('query', '').lower()
    selected_filters = request.args.getlist('filters')
    filters_str = ', '.join(selected_filters)
    limit = int(request.args.get('limit', 70))
    
    # For empty query and no filters, use the initial results endpoint
    if not query and not selected_filters:
        print("Empty search - redirecting to initial results")
        # Call the initial_results function directly but use the requested limit
        cached_result = search_cache.get('initial_results')
        if cached_result:
            print("Returning initial results from cache for empty search")
            return cached_result
            
        # If not cached, we'll just continue with the regular search flow
        # This will still be different from initial results due to the vector similarity
    
    # Create a cache key from the query parameters
    cache_key = f"{query}|{filters_str}|{limit}"
    cached_result = search_cache.get(cache_key)
    if cached_result:
        print(f"Returning cached search results for: {cache_key}")
        return cached_result
    
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
            abstract_text = row[10] if row[10] else "N/A"
            paper_data = {
                "Abstract": abstract_text,
                "description": abstract_text,  # Add description field for frontend compatibility
                "title": row[0] if row[0] else "Untitled Study",  # Add title field
                "Primary Outcome Area": row[7] if row[7] else "N/A",
                "Primary Outcome Measure": row[8] if row[8] else "N/A",
                "Treatment Duration": row[6] if row[6] else "N/A",
                "Publication Date": row[1].isoformat() if (row[1] and hasattr(row[1], 'isoformat')) else (str(row[1]) if row[1] else "N/A"),
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
            
        result = json.dumps(json_output, indent=4, sort_keys=False)
        # Cache the result before returning
        search_cache.set(cache_key, result)
        return result
        
    except Exception as e:
        print(f"Error executing vector search: {e}")
        return {}
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()

@app.route('/api/initial-results', methods=['GET'])
def get_initial_results():
    """Return initial set of results to preload on page load without waiting for user search"""
    print("Loading initial results")
    
    # Try to get from cache first
    cached_result = search_cache.get('initial_results')
    if cached_result:
        print("Returning initial results from cache")
        return cached_result
    
    limit = int(request.args.get('limit', 70))  # Use same default limit as search endpoint
    
    conn = get_db_connection()
    cursor = conn.cursor()
    
    # Use a more general query to fetch some initial data
    # This query just gets recent studies across treatments without any specific search terms
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
            0 AS distance,  -- No distance calculation needed for initial results
            spv.abstract
        FROM updated_treatment_data.semantic_paper_search_view spv
        JOIN updated_treatment_data.treatment_pubmed_link tpl ON spv.pmid = tpl.pmid
        JOIN updated_treatment_data.treatments t ON tpl.treatment_id = t.id
        ORDER BY spv.pub_date DESC  -- Get most recent studies
        LIMIT %s
    """
    
    try:
        cursor.execute(sql, [limit])
        results = cursor.fetchall()
        
        grouped_results = defaultdict(list)
        
        for row in results:
            abstract_text = row[10] if row[10] else "N/A"
            paper_data = {
                "Abstract": abstract_text,
                "description": abstract_text,  # Add description field for frontend compatibility
                "title": row[0] if row[0] else "Untitled Study",  # Add title field
                "Primary Outcome Area": row[7] if row[7] else "N/A",
                "Primary Outcome Measure": row[8] if row[8] else "N/A",
                "Treatment Duration": row[6] if row[6] else "N/A",
                "Publication Date": row[1].isoformat() if (row[1] and hasattr(row[1], 'isoformat')) else (str(row[1]) if row[1] else "N/A"),
                "Author": row[3] if row[3] else "N/A",
                "Study Title": row[0] if row[0] else "N/A",
                "PMID": str(row[2]) if row[2] else "N/A",
                "Full Text URL": row[4] if row[4] else None,
                "Similarity Score": 0.5,  # Default similarity score for initial results
                "Distance": 0.5  # Default distance for initial results
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
            
        # Also populate the medication filter while we're at it
        if len(available_filters["medication"]) == 0:
            # Extract unique treatment names and add them to the filter options
            available_treatments = list(set(treatment.lower() for treatment in grouped_results.keys() if treatment != "unknown"))
            available_filters["medication"] = sorted(available_treatments)
            print(f"Populated medication filters with {len(available_treatments)} treatments")
        
        result = json.dumps(json_output, indent=4, sort_keys=False)
        # Cache the result before returning
        search_cache.set('initial_results', result)
        return result
        
    except Exception as e:
        print(f"Error fetching initial results: {e}")
        return json.dumps([])
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()

if __name__ == '__main__':
    # Get port and debug settings from environment variables
    port = int(os.getenv('PYTHON_BACKEND_PORT', 5000))
    debug = os.getenv('DEBUG', 'true').lower() == 'true'
    app.run(debug=debug, port=port)