import pandas as pd

# File paths
files = [
    "Adjusted_ASD_Sheet_Info_V2.csv",
    "pubmed_papers_info.csv",
    "Individual_Treatment_Queries.csv"
]

# Load and combine
df_list = [pd.read_csv(f) for f in files]
combined = pd.concat(df_list, ignore_index=True)

# Drop duplicates by pmid
combined.drop_duplicates(subset='pmid', inplace=True)

# Extract pubmed_papers table
pubmed_columns = [
    'pmid', 'doi', 'title', 'pub_date', 'abstract',
    'authors', 'journal', 'keywords', 'url', 'affiliations'
]
pubmed_df = combined[pubmed_columns].copy()
pubmed_df.to_csv("cleaned_pubmed_papers.csv", index=False)

# Extract journals table
journals_df = combined[['journal', 'url']].drop_duplicates().copy()
journals_df['publisher'] = 'Unknown'
journals_df['impact_factor'] = None
journals_df['issn'] = None
journals_df = journals_df[['journal', 'publisher', 'impact_factor', 'issn', 'url']]
journals_df.to_csv("cleaned_journals.csv", index=False)

print("Exported cleaned_pubmed_papers.csv and cleaned_journals.csv")
