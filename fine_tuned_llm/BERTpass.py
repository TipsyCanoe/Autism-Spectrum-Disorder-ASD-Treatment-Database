from transformers import AutoModel, AutoTokenizer
from sentence_transformers import SentenceTransformer, models
import json
import re

input_path = "3data.json"




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

with open(input_path, 'r') as file:
    infile = json.load(file)

for i in range (201):
    embedding = sentence_model.encode(infile[i]['abstract'])
    embedding_list = embedding.tolist()
    infile[i]['vector'] = embedding_list

output_path = input_path.replace('.json', '_with_embeddings.json')
with open(output_path, 'w') as file:
    json.dump(infile, file, indent=4)

print(f"Embeddings added and saved to {output_path}")