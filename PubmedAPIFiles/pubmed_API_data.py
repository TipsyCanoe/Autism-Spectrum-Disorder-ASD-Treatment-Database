import pandas as pd
from pathlib import Path
from complete_query import importer

# Get terms from files
def get_terms(file_name, default):
    path = Path(file_name)
    terms = []
    if path.exists() and path.stat().st_size != 0:
        with open(file_name, "r") as file:
            for line in file:
                terms.append(line.strip())
    if not terms:
        terms = default
    return terms

def get_full_query():
    queries = []
    topic_queries = []
    full_query = []
    # Permanent keywords
    general_terms = ['Autism Spectrum Disorder[MeSH]', 'Autism[Title/Abstract]', 'Autistic Disorder[Title/Abstract]', 'ASD[Title/Abstract]']
    treatment_types_terms = ['randomized controlled trial[Publication Type]', 'randomized[Title/Abstract]', 'placebo[Title/Abstract]', '"Clinical Trial"[Publication Type]']

    # Default biologics and social terms
    biologics_terms = ['"Biological Products"[MeSH]', 'Biologic*[Title/Abstract]', '"Antibodies, monoclonal"[MeSH]', 'Cytokine Inhibitors[MeSH]', 'Immunotherapy[MeSH]']
    social_terms = ['rTMS[Title/Abstract]', '"Transcranial Magnetic Stimulation"[MeSH]', 'ECT[Title/Abstract]', '"Electroconvulsive Therapy"[MeSH]', 
    'Behavioral Therapy[MeSH]', '"Psychotherapy, Group"[MeSH]', 'Cognitive Therapy[MeSH]', 'Psychosocial Intervention[Title/Abstract]']

    biologics_terms = get_terms("biologic_terms.txt", biologics_terms)
    social_terms = get_terms("social_terms.txt", social_terms)

    # Building the query
    query_boxes = [general_terms, biologics_terms, social_terms, treatment_types_terms]

    for terms in query_boxes:
        topic_queries = ['{}'.format(topic) for topic in terms]
        queries.append('(' + ' OR '.join(topic_queries) + ')')

    full_query = queries[0] + ' AND (' + queries[1] + ' OR ' + queries[2] + ') AND ' + queries[3]
    return full_query

if __name__ == '__main__':
    full_query = get_full_query()

    # Importing the papers
    completePD = pd.DataFrame()
    completePD = importer.importPapers(completePD, full_query)
    completePD.to_excel('pubmed_papers_info.xlsx', index=False)