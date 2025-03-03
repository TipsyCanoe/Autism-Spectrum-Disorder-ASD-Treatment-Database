import faiss
import numpy as np
import sqlite3
from sentence_transformers import SentenceTransformer

### 1️⃣ Initialize the Embedding Model ###
model = SentenceTransformer('paraphrase-MiniLM-L6-v2')

# Example texts (replace with actual research papers)
texts = [
    "Artificial intelligence is transforming medicine.",
    "Deep learning outperforms traditional algorithms in image recognition.",
    "Neural networks achieve state-of-the-art results in NLP.",
    "AI is revolutionizing healthcare diagnostics.",
    "Natural language processing enables better communication between humans and machines."
]

### 2️⃣ Convert Text to Vectors ###
embeddings = model.encode(texts)
embeddings = np.array(embeddings).astype('float32')  # Ensure correct dtype for FAISS

### 3️⃣ Build FAISS Vector Index ###
dimension = embeddings.shape[1]  # Get vector size (e.g., 384)
index = faiss.IndexFlatL2(dimension)
index.add(embeddings)  # Add vectors to the FAISS index

### 4️⃣ Store Metadata (Original Strings) in SQLite ###
conn = sqlite3.connect('text_vector_db.db')
cursor = conn.cursor()

# Create a table for storing text and corresponding FAISS index
cursor.execute('''
    CREATE TABLE IF NOT EXISTS documents (
        id INTEGER PRIMARY KEY,
        text TEXT
    )
''')

# Insert text data with corresponding index
for i, text in enumerate(texts):
    cursor.execute('INSERT INTO documents (id, text) VALUES (?, ?)', (i, text))

conn.commit()

### 5️⃣ Perform a Similarity Search ###
query_text = "AI is transforming healthcare."
query_vector = model.encode([query_text]).astype('float32')

# Search FAISS for the most similar texts
distances, indices = index.search(query_vector, k=3)  # Find top 3 matches

# Fetch and display results from the SQLite database
print("\nQuery:", query_text)
print("Top Matches:")

for idx in indices[0]:
    cursor.execute('SELECT text FROM documents WHERE id = ?', (idx,))
    result = cursor.fetchone()

    if result:  # Ensure result is not None
        print(f"- {result[0]} (Index: {idx}, Distance: {distances[0][np.where(indices[0] == idx)][0]:.4f})")
    else:
        print(f"- No matching text found for Index: {idx}")

# Close database connection
conn.close()
