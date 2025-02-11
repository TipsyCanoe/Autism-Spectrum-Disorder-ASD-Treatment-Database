from flask import Flask, request, jsonify
from transformers import BertModel, BertTokenizer
import torch
import faiss
import numpy as np

app = Flask(__name__)

# Load BERT model and tokenizer
tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
model = BertModel.from_pretrained('bert-base-uncased')
model.eval()

# Load or create FAISS index
index = faiss.IndexFlatL2(768)  # 768 is the dimension of BERT embeddings

# Example data to add to the index (in practice, load your data)
example_texts = ["example sentence 1", "example sentence 2"]
example_embeddings = []

for text in example_texts:
    inputs = tokenizer(text, return_tensors='pt')
    with torch.no_grad():
        outputs = model(**inputs)
    embedding = outputs.last_hidden_state.mean(dim=1).numpy()
    example_embeddings.append(embedding)

example_embeddings = np.vstack(example_embeddings)
index.add(example_embeddings)

@app.route('/api/search', methods=['POST'])
def search():
    data = request.json
    selected_options = data['selectedOptions']
    
    # Process the selected options to create a query
    query = " ".join(selected_options)
    
    # Encode the query using BERT
    inputs = tokenizer(query, return_tensors='pt')
    with torch.no_grad():
        outputs = model(**inputs)
    query_embedding = outputs.last_hidden_state.mean(dim=1).numpy()
    
    # Search for nearest neighbors
    D, I = index.search(query_embedding, k=5)  # k is the number of nearest neighbors to return
    
    # In practice, map indices to actual data
    results = [example_texts[i] for i in I[0]]
    
    return jsonify(results)

if __name__ == '__main__':
    app.run(debug=True)