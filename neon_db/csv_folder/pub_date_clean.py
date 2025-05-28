import pandas as pd
from dateutil import parser

df = pd.read_csv("cleaned_pubmed_papers.csv")

def clean_date(val):
    try:
        return parser.parse(str(val), fuzzy=True).date()
    except Exception:
        return None

df['pub_date'] = df['pub_date'].apply(clean_date)
df.to_csv("cleaned_pubmed_papers_fixed.csv", index=False)
print("✅ Saved cleaned_pubmed_papers_fixed.csv with standardized dates.")
