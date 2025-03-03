from transformers import BertTokenizer, BertModel
import torch 
import matplotlib.pyplot as plt

# Load Med-BERT
MODEL_NAME = "bert-base-uncased"  # Adjust if Med-BERT has a specific name
tokenizer = BertTokenizer.from_pretrained(MODEL_NAME)
model = BertModel.from_pretrained(MODEL_NAME)

# Example medical text
text = "Patient diagnosed with anxiety and is a 12 year old female."

# Tokenize and convert to tensor
inputs = tokenizer(text, return_tensors="pt", truncation=True, padding=True, max_length=512)

# Generate embeddings
with torch.no_grad():
    outputs = model(**inputs)
    embeddings = outputs.last_hidden_state.mean(dim=1)  # Take mean for a single vector

print("Vector Shape:", embeddings.shape)  # Expected: (1, 768)
print("Vector Representation:", embeddings.numpy())

plt.hist(embeddings.numpy().flatten(), bins=50)
plt.xlabel("Embedding Value")
plt.ylabel("Frequency")
plt.title("Distribution of Embedding Values")
plt.show()