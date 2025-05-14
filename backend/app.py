from flask import Flask, request, jsonify
from flask_cors import CORS

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

resources = [
    {
        "id": 1,
        "title": "Understanding ASD in Children",
        "type": "Article",
        "publishDate": "2024-01-15",
        "author": "Dr. Jane Smith",
        "description": "A comprehensive overview of autism spectrum disorder in children and best practices for support.",
        "tags": ["autism", "children", "support"],
        "url": "https://example.com/article1",
        "age": "6-12",
        "symptom": "social",
        "gender": "male"
    },
    {
        "id": 2,
        "title": "Managing ADHD Symptoms in Adolescents",
        "type": "Guide",
        "publishDate": "2023-11-20",
        "author": "Dr. Michael Johnson",
        "description": "Strategies for parents and educators to help adolescents manage ADHD symptoms.",
        "tags": ["adhd", "adolescents", "management"],
        "url": "https://example.com/guide1",
        "age": "13-17",
        "symptom": "adhd",
        "gender": "female"
    },
    {
        "id": 3,
        "title": "Early Intervention Techniques for ASD",
        "type": "Research Paper",
        "publishDate": "2024-02-05",
        "author": "Dr. Sarah Williams",
        "description": "Research findings on effective early intervention techniques for children with autism spectrum disorder.",
        "tags": ["autism", "early intervention", "research"],
        "url": "https://example.com/research1",
        "age": "0-5",
        "symptom": "social",
        "gender": "nonbinary"
    },
    {
        "id": 4,
        "title": "Understanding Hyperactivity in Young Adults",
        "type": "Article",
        "publishDate": "2023-12-10",
        "author": "Dr. Robert Chen",
        "description": "An exploration of hyperactivity in young adults and its impact on daily functioning.",
        "tags": ["hyperactivity", "young adults"],
        "url": "https://example.com/article2",
        "age": "18-25",
        "symptom": "hyperactivity",
        "gender": "male"
    },
    {
        "id": 5,
        "title": "Anxiety Management Strategies for Children with ASD",
        "type": "Guide",
        "publishDate": "2024-03-01",
        "author": "Dr. Emily Davis",
        "description": "Practical strategies to help children with ASD manage anxiety and reactivity.",
        "tags": ["anxiety", "children", "autism"],
        "url": "https://example.com/guide2",
        "age": "6-12",
        "symptom": "anxiety-reactivity",
        "gender": "female"
    }
]


available_filters = {
    "age": ["0-5", "6-12", "13-17", "18-25", "26-64", "65+"],
    "symptom": [
        "irritability", "adhd", "hyperactivity", "social", 
        "attention-hyperactivity", "asd-severity", 
        "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance", 
        "anxiety-reactivity"
    ],
    "gender": ["male", "female", "nonbinary"]
}

@app.route('/api/filters', methods=['GET'])
def get_filters():
    """Return available filter options"""
    return jsonify(available_filters)

@app.route('/api/search', methods=['GET'])
def search():
    """Search resources based on query and filters"""
    query = request.args.get('query', '').lower()
    embedding = sentence_model.encode(query)
    print(embedding)
    print(query)
    selected_filters = request.args.getlist('filters')
    parsed_filters = {}
    for filter_str in selected_filters:
        if ':' in filter_str:
            category, value = filter_str.split(':', 1)
            if category not in parsed_filters:
                parsed_filters[category] = []
            parsed_filters[category].append(value)
    
    
    filtered_results = resources
    
    
    for category, values in parsed_filters.items():
        if values:
            filtered_results = [
                resource for resource in filtered_results
                if resource.get(category) in values
            ]
    
    return jsonify({"results": filtered_results, "total": len(filtered_results)})
if __name__ == '__main__':
    app.run(debug=True, port=5000)