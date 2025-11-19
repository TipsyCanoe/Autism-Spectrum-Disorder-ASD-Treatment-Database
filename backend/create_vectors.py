import pandas as pd
import csv
from transformers import AutoTokenizer 
from sentence_transformers import SentenceTransformer, models
from tqdm import tqdm

input_csv = "FullLLMResults.csv"
output_csv = "LLM_with_embeddings.csv"
model_name = "Charangan/MedBERT"
df = pd.read_csv(input_csv, encoding="utf-8", engine="python")

if "abstract" not in df.columns:
    raise ValueError("CSV must contain an 'abstract' column.")

tokenizer = AutoTokenizer.from_pretrained(model_name)
word_embedding_model = models.Transformer(model_name)
pooling_model = models.Pooling(
    word_embedding_model.get_word_embedding_dimension(),
    pooling_mode_mean_tokens=True,
    pooling_mode_cls_token=False,
    pooling_mode_max_tokens=False
)
sentence_model = SentenceTransformer(modules=[word_embedding_model, pooling_model])

embeddings = []
for text in tqdm(df["abstract"].fillna(""), desc="Encoding abstracts"):
    if isinstance(text, str) and text.strip():
        emb = sentence_model.encode(text)
        embeddings.append(emb.tolist())
    else:
        embeddings.append(None)

df["vector"] = embeddings
df["AI"] = True

df["vector"] = df["vector"].apply(lambda v: str(v) if v is not None else "")

df.to_csv(output_csv, index=False, encoding="utf-8")
