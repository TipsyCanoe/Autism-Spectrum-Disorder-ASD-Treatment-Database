import pandas as pd
from pathlib import Path
from complete_query import importer
import re
import os
from datetime import datetime

def get_treatments_queries(treatments):
    treatment_queries = []
    if treatments:
        for i in range(len(treatments)):
            for keyword in ["AND", "OR"]:
                treatments[i] = re.sub(fr"\b{keyword}\b", keyword, treatments[i], flags=re.IGNORECASE)
        treatment_queries = ['({})'.format(treatment) for treatment in treatments]
    return treatment_queries

def get_full_treatment_queries(treatment_queries, queries):
    full_queries = []
    for treatment in treatment_queries:
        full_queries.append(treatment + ' AND ' + queries[0] + ' AND ' + queries[1])
    return full_queries

if __name__ == '__main__':
    # Use path relative to this script
    path = str(Path(__file__).parent / 'Treatment_Names.xlsx')

    # Null check
    if not os.path.isfile(path):
        print(f"File not found: {path}; skipping")
        exit(1)

    oldDf = pd.read_excel(path, 'Sheet1', usecols="A")

    treatments = oldDf['Treatment name'].tolist() 

    queries = []

    full_queries = []

    treatment_queries = get_treatments_queries(treatments)

    general_terms = ['Autism Spectrum Disorder[MeSH]', 'Autism[Title/Abstract]', 'Autistic Disorder[Title/Abstract]', 'ASD[Title/Abstract]']

    treatment_types_terms = ['"randomized controlled trial"[Publication Type]', 'randomized controlled trial[Title/Abstract]']

    for terms in [general_terms, treatment_types_terms]:
        topic_queries = ['{}'.format(topic) for topic in terms]
        queries.append('(' + ' OR '.join(topic_queries) + ')')
    
    full_queries = get_full_treatment_queries(treatment_queries, queries)

    completePD = pd.DataFrame()

    for query in full_queries:
        completePD = importer.importPapers(completePD, query)

    today = datetime.now().strftime("%Y-%m-%d")
    completePD.to_excel(f"pubmed_treatment_info_{today}.xlsx", index=False)