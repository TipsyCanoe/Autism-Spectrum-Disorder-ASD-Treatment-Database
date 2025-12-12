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


embedding = sentence_model.encode("patient reports shortness of breath and chest pain.")

#print(embedding)

#from sklearn.metrics.pairwise import cosine_similarity

#this prints a similarity matrix the last one should be the least similar
"""
sentences = [
    "The patient has a fever and cough.",
    "The individual is experiencing a high temperature and a sore throat.",
    "I like playing video games.",
]

embeddings = sentence_model.encode(sentences)
similarity_matrix = cosine_similarity(embeddings)

print(similarity_matrix)
"""

#this is to test which is closest to the search query should be the first
"""
documents = [
    "Chest X-ray shows signs of pneumonia.",
    "The patient was prescribed antibiotics.",
    "I like pasta.",
]

query = "Lung infection visible in imaging"

doc_embeddings = sentence_model.encode(documents)
query_embedding = sentence_model.encode(query)

# Get cosine similarities
sims = cosine_similarity([query_embedding], doc_embeddings)
print(sims)
"""