from flask import Flask, request, jsonify
from flask_cors import CORS
import re

app = Flask(__name__)
CORS(app) 

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

resources = {
    "aripiprazole": [
        {
        "Study Title": "Efficacy of Aripiprazole in Adolescents with Bipolar Disorder",
        "Duration": "8 weeks",
        "n": 100,
        "M:F ratio": "60:40",
        "Age range/mean": "13-17 / 15",
        "Medication/Treatment Dose Range": "2–30 mg/day",
        "Primary Outcome Area": "Mania symptoms",
        "Primary Outcome Measures": "YMRS",
        "Results: Primary measure": "Significant reduction in YMRS scores",
        "Secondary Outcome Area": "Depressive symptoms",
        "Secondary Outcome Measures": "CDRS-R",
        "Results: Secondary Measures": "Mild improvement",
        "Tolerability/Side Effects": "Akathisia, nausea",
        "Safety": "No severe adverse events",
        "Drop Out Rate": "15%"
        },
        {
        "Study Title": "Long-Term Safety of Aripiprazole in Pediatric Schizophrenia",
        "Duration": "24 weeks",
        "n": 150,
        "M:F ratio": "70:80",
        "Age range/mean": "10-17 / 14.2",
        "Medication/Treatment Dose Range": "5–20 mg/day",
        "Primary Outcome Area": "Psychotic symptoms",
        "Primary Outcome Measures": "PANSS",
        "Results: Primary measure": "Reduction in PANSS positive subscale",
        "Secondary Outcome Area": "Cognitive function",
        "Secondary Outcome Measures": "CANTAB",
        "Results: Secondary Measures": "No significant change",
        "Tolerability/Side Effects": "Weight gain, fatigue",
        "Safety": "No serious adverse events reported",
        "Drop Out Rate": "12%"
        },
        {
        "Study Title": "Aripiprazole Augmentation in Adolescent Depression",
        "Duration": "10 weeks",
        "n": 80,
        "M:F ratio": "30:50",
        "Age range/mean": "14-18 / 16",
        "Medication/Treatment Dose Range": "2–10 mg/day",
        "Primary Outcome Area": "Treatment-resistant depression",
        "Primary Outcome Measures": "HAM-D, CGI-I",
        "Results: Primary measure": "Improved response rate vs. placebo",
        "Secondary Outcome Area": "Functioning",
        "Secondary Outcome Measures": "GAF",
        "Results: Secondary Measures": "Moderate improvement",
        "Tolerability/Side Effects": "Insomnia, irritability",
        "Safety": "Well-monitored; no significant safety concerns",
        "Drop Out Rate": "18%"
        }
    ],
    "citalopram": [
        {
        "Study Title": "Citalopram for Major Depressive Disorder in Adolescents",
        "Duration": "12 weeks",
        "n": 120,
        "M:F ratio": "50:70",
        "Age range/mean": "12-18 / 15.3",
        "Medication/Treatment Dose Range": "10–40 mg/day",
        "Primary Outcome Area": "Depression",
        "Primary Outcome Measures": "HAM-D, CGI",
        "Results: Primary measure": "Marked reduction in HAM-D scores",
        "Secondary Outcome Area": "Anxiety symptoms",
        "Secondary Outcome Measures": "SCARED",
        "Results: Secondary Measures": "Moderate improvement",
        "Tolerability/Side Effects": "GI upset, headache",
        "Safety": "Well-tolerated",
        "Drop Out Rate": "10%"
        },
        {
        "Study Title": "Effectiveness of Citalopram in Adolescent Anxiety Disorders",
        "Duration": "10 weeks",
        "n": 90,
        "M:F ratio": "40:50",
        "Age range/mean": "13-17 / 15.1",
        "Medication/Treatment Dose Range": "10–20 mg/day",
        "Primary Outcome Area": "Generalized anxiety",
        "Primary Outcome Measures": "PARS",
        "Results: Primary measure": "Significant reduction in anxiety symptoms",
        "Secondary Outcome Area": "School performance",
        "Secondary Outcome Measures": "Teacher ratings",
        "Results: Secondary Measures": "Slight improvement",
        "Tolerability/Side Effects": "Somnolence, dry mouth",
        "Safety": "No serious concerns",
        "Drop Out Rate": "8%"
        },
        {
        "Study Title": "Citalopram in Youth with Autism and Comorbid Depression",
        "Duration": "16 weeks",
        "n": 60,
        "M:F ratio": "45:15",
        "Age range/mean": "10-15 / 12.5",
        "Medication/Treatment Dose Range": "5–20 mg/day",
        "Primary Outcome Area": "Depressive symptoms",
        "Primary Outcome Measures": "CDRS-R",
        "Results: Primary measure": "Improvement in mood symptoms",
        "Secondary Outcome Area": "Repetitive behaviors",
        "Secondary Outcome Measures": "RBS-R",
        "Results: Secondary Measures": "No significant change",
        "Tolerability/Side Effects": "Agitation, appetite changes",
        "Safety": "Well tolerated with close monitoring",
        "Drop Out Rate": "13%"
        }
    ],
    "fluoxetine": [
        {
        "Study Title": "Fluoxetine in the Treatment of Pediatric Major Depression",
        "Duration": "12 weeks",
        "n": 150,
        "M:F ratio": "70:80",
        "Age range/mean": "8-18 / 14",
        "Medication/Treatment Dose Range": "10–40 mg/day",
        "Primary Outcome Area": "Depressive symptoms",
        "Primary Outcome Measures": "CDRS-R, CGI",
        "Results: Primary measure": "Significant symptom reduction compared to placebo",
        "Secondary Outcome Area": "Functioning",
        "Secondary Outcome Measures": "Children’s Global Assessment Scale (CGAS)",
        "Results: Secondary Measures": "Improved functioning",
        "Tolerability/Side Effects": "Headache, sleep disturbances",
        "Safety": "No suicidality increase noted",
        "Drop Out Rate": "11%"
        },
        {
        "Study Title": "Long-Term Outcomes of Fluoxetine in Adolescents with Depression",
        "Duration": "6 months",
        "n": 100,
        "M:F ratio": "45:55",
        "Age range/mean": "12-17 / 15.2",
        "Medication/Treatment Dose Range": "20–60 mg/day",
        "Primary Outcome Area": "Depression relapse prevention",
        "Primary Outcome Measures": "HAM-D",
        "Results: Primary measure": "Lower relapse rate vs. placebo",
        "Secondary Outcome Area": "Academic performance",
        "Secondary Outcome Measures": "Teacher reports",
        "Results: Secondary Measures": "No significant difference",
        "Tolerability/Side Effects": "Dry mouth, mild agitation",
        "Safety": "Stable vitals and labs",
        "Drop Out Rate": "9%"
        },
        {
        "Study Title": "Fluoxetine for Adolescent Obsessive-Compulsive Disorder",
        "Duration": "14 weeks",
        "n": 85,
        "M:F ratio": "40:45",
        "Age range/mean": "13-17 / 14.8",
        "Medication/Treatment Dose Range": "10–60 mg/day",
        "Primary Outcome Area": "OCD symptoms",
        "Primary Outcome Measures": "CY-BOCS",
        "Results: Primary measure": "Marked improvement in CY-BOCS scores",
        "Secondary Outcome Area": "Anxiety",
        "Secondary Outcome Measures": "SCARED",
        "Results: Secondary Measures": "Moderate reduction in anxiety",
        "Tolerability/Side Effects": "GI discomfort, jitteriness",
        "Safety": "No major safety events",
        "Drop Out Rate": "12%"
        }
    ],
    "risperidone": [
        {
        "Study Title": "Risperidone in Adolescents with Schizophrenia",
        "Duration": "8 weeks",
        "n": 110,
        "M:F ratio": "65:45",
        "Age range/mean": "13-18 / 15.5",
        "Medication/Treatment Dose Range": "0.5–6 mg/day",
        "Primary Outcome Area": "Psychotic symptoms",
        "Primary Outcome Measures": "PANSS",
        "Results: Primary measure": "Significant reduction in total PANSS score",
        "Secondary Outcome Area": "Negative symptoms",
        "Secondary Outcome Measures": "PANSS negative subscale",
        "Results: Secondary Measures": "Moderate improvement",
        "Tolerability/Side Effects": "Weight gain, sedation",
        "Safety": "Prolactin elevation observed",
        "Drop Out Rate": "14%"
        },
        {
        "Study Title": "Risperidone for Irritability in Autism Spectrum Disorders",
        "Duration": "10 weeks",
        "n": 96,
        "M:F ratio": "75:21",
        "Age range/mean": "5-16 / 9.7",
        "Medication/Treatment Dose Range": "0.25–3.5 mg/day",
        "Primary Outcome Area": "Irritability and aggression",
        "Primary Outcome Measures": "ABC Irritability subscale",
        "Results: Primary measure": "Substantial improvement vs. placebo",
        "Secondary Outcome Area": "Repetitive behaviors",
        "Secondary Outcome Measures": "RBS-R",
        "Results: Secondary Measures": "Slight improvement",
        "Tolerability/Side Effects": "Drowsiness, increased appetite",
        "Safety": "Some metabolic changes monitored",
        "Drop Out Rate": "7%"
        },
        {
        "Study Title": "Risperidone for Bipolar Mania in Youth",
        "Duration": "6 weeks",
        "n": 130,
        "M:F ratio": "68:62",
        "Age range/mean": "10-17 / 14.6",
        "Medication/Treatment Dose Range": "0.5–4 mg/day",
        "Primary Outcome Area": "Manic symptoms",
        "Primary Outcome Measures": "YMRS",
        "Results: Primary measure": "Significant symptom control",
        "Secondary Outcome Area": "Global functioning",
        "Secondary Outcome Measures": "CGAS",
        "Results: Secondary Measures": "Improved functioning scores",
        "Tolerability/Side Effects": "Mild EPS, somnolence",
        "Safety": "No serious adverse effects",
        "Drop Out Rate": "13%"
        }
    ],
    "sertraline": [
        {
        "Study Title": "Sertraline for Pediatric Anxiety Disorders",
        "Duration": "12 weeks",
        "n": 145,
        "M:F ratio": "70:75",
        "Age range/mean": "7-17 / 12.8",
        "Medication/Treatment Dose Range": "25–200 mg/day",
        "Primary Outcome Area": "Generalized anxiety, separation anxiety, social phobia",
        "Primary Outcome Measures": "PARS",
        "Results: Primary measure": "Significant reduction in anxiety scores",
        "Secondary Outcome Area": "Parent-rated outcomes",
        "Secondary Outcome Measures": "Child Behavior Checklist (CBCL)",
        "Results: Secondary Measures": "Improved behavior ratings",
        "Tolerability/Side Effects": "GI issues, insomnia",
        "Safety": "No suicidality noted",
        "Drop Out Rate": "9%"
        },
        {
        "Study Title": "Sertraline for OCD in Adolescents",
        "Duration": "10 weeks",
        "n": 90,
        "M:F ratio": "40:50",
        "Age range/mean": "13-18 / 15",
        "Medication/Treatment Dose Range": "50–200 mg/day",
        "Primary Outcome Area": "OCD symptoms",
        "Primary Outcome Measures": "CY-BOCS",
        "Results: Primary measure": "Moderate to strong improvement",
        "Secondary Outcome Area": "Family burden",
        "Secondary Outcome Measures": "Family Accommodation Scale",
        "Results: Secondary Measures": "Slight improvement",
        "Tolerability/Side Effects": "Headache, activation",
        "Safety": "Well-tolerated",
        "Drop Out Rate": "10%"
        },
        {
        "Study Title": "Sertraline for Depression in Adolescents Post-Trauma",
        "Duration": "8 weeks",
        "n": 75,
        "M:F ratio": "35:40",
        "Age range/mean": "12-16 / 14",
        "Medication/Treatment Dose Range": "25–100 mg/day",
        "Primary Outcome Area": "Depression",
        "Primary Outcome Measures": "CDRS-R",
        "Results: Primary measure": "Symptom reduction, especially in moderate-to-severe cases",
        "Secondary Outcome Area": "PTSD symptoms",
        "Secondary Outcome Measures": "Child PTSD Symptom Scale",
        "Results: Secondary Measures": "Mild reduction",
        "Tolerability/Side Effects": "Fatigue, occasional nausea",
        "Safety": "Safe with monitoring",
        "Drop Out Rate": "12%"
        }
    ]
}

available_filters = {
    "age": ["0-5", "6-12", "13-17", "18-25", "26-64", "65+"],
    "symptom": [
        "irritability", "adhd", "hyperactivity", "social", 
        "attention-hyperactivity", "asd-severity", 
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
            "asd-severity": ["autism", "asd"],
            "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance": ["lethargy", "withdrawal", "stereotypy"],
            "anxiety-reactivity": ["anxiety", "anxious", "reactivity"]
        }
        
        terms = symptom_terms.get(symptom, [symptom])
        if any(term in searchable_text for term in terms):
            return True
    
    return False

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
    """Search resources based on query and filters"""
    query = request.args.get('query', '').strip()
    selected_filters = request.args.getlist('filters')
    
    # Parse the filters
    parsed_filters = {}
    for filter_str in selected_filters:
        if ':' in filter_str:
            category, value = filter_str.split(':', 1)
            if category not in parsed_filters:
                parsed_filters[category] = []
            parsed_filters[category].append(value)
    
    print(f"Query: {query}")
    print(f"Parsed filters: {parsed_filters}")
    
    # If query is provided, generate embedding for semantic search
    if query:
        embedding = sentence_model.encode(query)
        print(f"Generated embedding shape: {embedding.shape}")
    
    # Start with all resources
    filtered_results = {}
    
    medications_to_search = parsed_filters.get('medication', list(resources.keys()))
    
    for medication in medications_to_search:
        if medication not in resources:
            continue
            
        filtered_studies = []
        
        for study in resources[medication]:
            # Check all filter conditions
            passes_filters = True
            
            if 'age' in parsed_filters:
                if not matches_age_filter(study.get('Age range/mean'), parsed_filters['age']):
                    passes_filters = False
            
            if 'gender' in parsed_filters:
                if not matches_gender_filter(study.get('M:F ratio'), parsed_filters['gender']):
                    passes_filters = False
            
            if 'symptom' in parsed_filters:
                if not matches_symptom_filter(study, parsed_filters['symptom']):
                    passes_filters = False
            
            
            if passes_filters:
                filtered_studies.append(study)
        
        # Only include medication if it has matching studies
        if filtered_studies:
            filtered_results[medication] = filtered_studies
    
    # Count total studies
    total_studies = sum(len(studies) for studies in filtered_results.values())
    
    print(f"Returning {total_studies} studies across {len(filtered_results)} medications")
    
    return jsonify({
        "results": filtered_results, 
        "total": total_studies
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)