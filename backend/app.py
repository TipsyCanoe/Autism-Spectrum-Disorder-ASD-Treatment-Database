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
    "symptom": [],  # Will be populated from primary_outcome_area column
    "gender": ["male", "female", "nonbinary"],
    "medication": []  # Will be populated from treatment_name column
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
    connection_url = os.getenv('DATABASE_URL', NEON_DATABASE_URL)
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
    
    # Load medication and symptom options from database if they're empty
    if not available_filters["medication"] or not available_filters["symptom"]:
        try:
            conn = get_db_connection()
            cursor = conn.cursor()
            
            # Get distinct treatment names for medication filter
            if not available_filters["medication"]:
                cursor.execute("""
                    SELECT DISTINCT treatment_name 
                    FROM jim_data.data_embedded 
                    WHERE treatment_name IS NOT NULL AND treatment_name <> '' 
                    ORDER BY treatment_name ASC
                """)
                medications = [row[0].lower() for row in cursor.fetchall()]
                available_filters["medication"] = medications
                print(f"Loaded {len(medications)} medications for filter")
            
            # Get distinct primary outcome areas for symptom filter
            if not available_filters["symptom"]:
                cursor.execute("""
                    SELECT DISTINCT primary_outcome_area 
                    FROM jim_data.data_embedded 
                    WHERE primary_outcome_area IS NOT NULL AND primary_outcome_area <> '' 
                    ORDER BY primary_outcome_area ASC
                """)
                symptoms = [row[0] for row in cursor.fetchall()]
                available_filters["symptom"] = symptoms
                print(f"Loaded {len(symptoms)} primary outcome areas for symptom filter")
                
        except Exception as e:
            print(f"Error loading filters from database: {e}")
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
    limit = int(request.args.get('limit', 10000))  # Increased default limit
    
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
            title,
            pub_date,
            pmid,
            authors,
            url,
            treatment_name,
            study_duration,
            primary_outcome_area,
            primary_outcome_measures,
            vector::vector <=> %s::vector AS distance,
            abstract,
            authors,
            pub_date,
            study_type,
            sample_size,
            male_female_ratio,
            age_range,
            treatment_dose_range,
            results_primary_measure,
            secondary_outcome_area,
            secondary_outcome_measures,
            tolerability_side_effects,
            safety,
            drop_out_rate,
            ethnicity_percentages,
            NULL as notes,
            NULL as sequence_generation_selection_bias,
            NULL as allocation_concealment_selection_bias,
            NULL as outcome_assessors_blinding_detection_bias,
            NULL as clinician_and_participant_blinding_performance_bias,
            NULL as incomplete_outcome_data_attrition_bias,
            NULL as selective_outcome_reporting_reporting_bias,
            NULL as notes_on_biases,
            ai,
            NULL as age_min,
            NULL as age_max,
            NULL as males_in_study,
            NULL as females_in_study,
            journal,
            affiliations
        FROM jim_data.data_embedded
        WHERE vector IS NOT NULL AND vector::vector <=> %s::vector < %s
        ORDER BY distance ASC
        LIMIT %s
    """
    
    params = [query_embedding.tolist(), query_embedding.tolist(), similarity_threshold, limit]
    
    try:
        cursor.execute(sql, params)
        results = cursor.fetchall()
        
        grouped_results = defaultdict(list)
        
        for row in results:
            abstract_text = row[10] if row[10] else "Not specified in article"
            paper_data = {
                "Abstract": abstract_text,
                "description": abstract_text,  # Add description field for frontend compatibility
                "title": row[0] if row[0] else "Untitled Study",
                "Primary Outcome Area": row[7] if row[7] else "Not specified in article",
                "Primary Outcome Measure": row[8] if row[8] else "Not specified in article",
                "Treatment Duration": row[6] if row[6] else "Not specified in article",
                "Publication Date": row[1].isoformat() if (row[1] and hasattr(row[1], 'isoformat')) else (str(row[1]) if row[1] else "Not specified in article"),
                "Author": row[11] if row[11] else (row[3] if row[3] else "Not specified in article"),
                "Study Title": row[0] if row[0] else "Not specified in article",
                "PMID": str(row[2]) if row[2] else "Not specified in article",
                "Full Text URL": row[4] if row[4] else None,
                "Similarity Score": round(1 - row[9], 3) if row[9] is not None else 0,
                "Distance": round(row[9], 3) if row[9] is not None else 1,
                # Additional fields from new schema
                "Study Type": row[13] if row[13] else "Not specified in article",
                "Sample Size": row[14] if row[14] else "Not specified in article",
                "M:F Ratio": row[15] if row[15] else "Not specified in article",
                "Age Range (years)": row[16] if row[16] else "Not specified in article",
                "Medication/Treatment Dose Range": row[17] if row[17] else "Not specified in article",
                "Results: Primary measure": row[18] if row[18] else "Not specified in article",
                "Secondary Outcome Area": row[19] if row[19] else "Not specified in article",
                "Secondary Outcome Measures": row[20] if row[20] else "Not specified in article",
                "Tolerability/Side Effects": row[21] if row[21] else "Not specified in article",
                "Safety": row[22] if row[22] else "Not specified in article",
                "Drop Out Rate": row[23] if row[23] else "Not specified in article",
                "Race/Ethnicity Percentages": row[24] if row[24] else "Not specified in article",
                "Notes": row[25] if row[25] else "Not specified in article",
                "Sequence Generation (selection bias)": row[26] if row[26] else "Not specified in article",
                "Allocation Concealment (selection bias)": row[27] if row[27] else "Not specified in article",
                "Outcome Assessors Blinding (detection bias)": row[28] if row[28] else "Not specified in article",
                "Clinician and Participant Blinding (performance bias)": row[29] if row[29] else "Not specified in article",
                "Incomplete outcome data (attrition bias)": row[30] if row[30] else "Not specified in article",
                "Selective outcome reporting (reporting bias)": row[31] if row[31] else "Not specified in article",
                "Notes on Biases": row[32] if row[32] else "Not specified in article",
                "AI": row[33] if row[33] else "Not specified in article",
                "age_min": row[34] if row[34] is not None else "Not specified in article",
                "age_max": row[35] if row[35] is not None else "Not specified in article",
                "males_in_study": row[36] if row[36] is not None else "Not specified in article",
                "females_in_study": row[37] if row[37] is not None else "Not specified in article",
                "Journal": row[38] if row[38] else "Not specified in article",
                "Commercial Affiliation": row[39] if row[39] else "Not specified in article"
            }
            
            treatment_name = row[5].lower() if row[5] else "unknown"
            grouped_results[treatment_name].append(paper_data)
        
        # Sort studies within each treatment by similarity score (best first)
        for treatment_name in grouped_results:
            grouped_results[treatment_name].sort(key=lambda x: x["Similarity Score"], reverse=True)
        
        # Sort treatments by best similarity score (treatment with highest scoring study first)
        final_results = dict(grouped_results)
        final_results = dict(sorted(final_results.items(), 
                                  key=lambda x: max([s["Similarity Score"] for s in x[1]], default=0), 
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
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
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
    
    limit = int(request.args.get('limit', 10000))  # Increased default limit
    
    conn = get_db_connection()
    cursor = conn.cursor()
    
    # Use a more general query to fetch some initial data
    # This query just gets recent studies across treatments without any specific search terms
    sql = """
        SELECT 
            title,
            pub_date,
            pmid,
            authors,
            url,
            treatment_name,
            study_duration,
            primary_outcome_area,
            primary_outcome_measures,
            0 AS distance,  -- No distance calculation needed for initial results
            abstract,
            authors,
            pub_date,
            study_type,
            sample_size,
            male_female_ratio,
            age_range,
            treatment_dose_range,
            results_primary_measure,
            secondary_outcome_area,
            secondary_outcome_measures,
            tolerability_side_effects,
            safety,
            drop_out_rate,
            ethnicity_percentages,
            NULL as notes,
            NULL as sequence_generation_selection_bias,
            NULL as allocation_concealment_selection_bias,
            NULL as outcome_assessors_blinding_detection_bias,
            NULL as clinician_and_participant_blinding_performance_bias,
            NULL as incomplete_outcome_data_attrition_bias,
            NULL as selective_outcome_reporting_reporting_bias,
            NULL as notes_on_biases,
            ai,
            NULL as age_min,
            NULL as age_max,
            NULL as males_in_study,
            NULL as females_in_study,
            journal,
            affiliations
        FROM jim_data.data_embedded
        ORDER BY pub_date DESC NULLS LAST  -- Get most recent studies
        LIMIT %s
    """
    
    try:
        cursor.execute(sql, [limit])
        results = cursor.fetchall()
        
        grouped_results = defaultdict(list)
        
        for row in results:
            abstract_text = row[10] if row[10] else "Not specified in article"
            paper_data = {
                "Abstract": abstract_text,
                "description": abstract_text,  # Add description field for frontend compatibility
                "title": row[0] if row[0] else "Untitled Study",
                "Primary Outcome Area": row[7] if row[7] else "Not specified in article",
                "Primary Outcome Measure": row[8] if row[8] else "Not specified in article",
                "Treatment Duration": row[6] if row[6] else "Not specified in article",
                "Publication Date": row[1].isoformat() if (row[1] and hasattr(row[1], 'isoformat')) else (str(row[1]) if row[1] else "Not specified in article"),
                "Author": row[11] if row[11] else (row[3] if row[3] else "Not specified in article"),
                "Study Title": row[0] if row[0] else "Not specified in article",
                "PMID": str(row[2]) if row[2] else "Not specified in article",
                "Full Text URL": row[4] if row[4] else None,
                "Similarity Score": 0.5,  # Default similarity score for initial results
                "Distance": 0.5,  # Default distance for initial results
                # Additional fields from new schema
                "Study Type": row[13] if row[13] else "Not specified in article",
                "Sample Size": row[14] if row[14] else "Not specified in article",
                "M:F Ratio": row[15] if row[15] else "Not specified in article",
                "Age Range (years)": row[16] if row[16] else "Not specified in article",
                "Medication/Treatment Dose Range": row[17] if row[17] else "Not specified in article",
                "Results: Primary measure": row[18] if row[18] else "Not specified in article",
                "Secondary Outcome Area": row[19] if row[19] else "Not specified in article",
                "Secondary Outcome Measures": row[20] if row[20] else "Not specified in article",
                "Tolerability/Side Effects": row[21] if row[21] else "Not specified in article",
                "Safety": row[22] if row[22] else "Not specified in article",
                "Drop Out Rate": row[23] if row[23] else "Not specified in article",
                "Race/Ethnicity Percentages": row[24] if row[24] else "Not specified in article",
                "Notes": row[25] if row[25] else "Not specified in article",
                "Sequence Generation (selection bias)": row[26] if row[26] else "Not specified in article",
                "Allocation Concealment (selection bias)": row[27] if row[27] else "Not specified in article",
                "Outcome Assessors Blinding (detection bias)": row[28] if row[28] else "Not specified in article",
                "Clinician and Participant Blinding (performance bias)": row[29] if row[29] else "Not specified in article",
                "Incomplete outcome data (attrition bias)": row[30] if row[30] else "Not specified in article",
                "Selective outcome reporting (reporting bias)": row[31] if row[31] else "Not specified in article",
                "Notes on Biases": row[32] if row[32] else "Not specified in article",
                "AI": row[33] if row[33] else "Not specified in article",
                "age_min": row[34] if row[34] is not None else "Not specified in article",
                "age_max": row[35] if row[35] is not None else "Not specified in article",
                "males_in_study": row[36] if row[36] is not None else "Not specified in article",
                "females_in_study": row[37] if row[37] is not None else "Not specified in article",
                "Journal": row[38] if row[38] else "Not specified in article",
                "Commercial Affiliation": row[39] if row[39] else "Not specified in article"
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

@app.route('/api/clear-cache', methods=['POST'])
def clear_cache():
    """Clear all caches"""
    search_cache.cache.clear()
    filter_cache.cache.clear()
    return jsonify({"message": "Cache cleared successfully"}), 200

if __name__ == '__main__':
    # Get port and debug settings from environment variables
    port = int(os.getenv('PYTHON_BACKEND_PORT', 5000))
    debug = os.getenv('DEBUG', 'true').lower() == 'true'
    app.run(debug=debug, port=port)