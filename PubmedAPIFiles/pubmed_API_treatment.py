import pandas as pd
from pathlib import Path
from complete_query import importer

# May need to adjust path depending on what directory you run this in
path = str(Path.cwd()) + '/../PubmedAPIFiles/Treatment_Names.xlsx'

oldDf = pd.read_excel(path, 'Sheet1', usecols="A")

treatments = oldDf['Treatment name'].tolist() 

queries = []
treatment_queries = []
full_queries = []

if treatments:
    for i in range(len(treatments)):
        if ' and ' in treatments[i] or ' or ' in treatments[i]:
            treatments[i] = treatments[i].replace(' and ', ' AND ')
            treatments[i] = treatments[i].replace(' or ', ' OR ') 
    treatment_queries = ['({})'.format(treatment) for treatment in treatments]

general_terms = ['Autism Spectrum Disorder[MeSH]', 'Autism[Title/Abstract]', 'Autistic Disorder[Title/Abstract]', 'ASD[Title/Abstract]']

treatment_types_terms = ['"randomized controlled trial"[Publication Type]', 'randomized controlled trial[Title/Abstract]']

for terms in [general_terms, treatment_types_terms]:
    topic_queries = ['{}'.format(topic) for topic in terms]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

for treatment in treatment_queries:
    full_queries.append(treatment + ' AND ' + queries[0] + ' AND ' + queries[1])
# print(full_queries[0])

completePD = pd.DataFrame()

count = 0

for query in full_queries:
    completePD = importer.importPapers(completePD, query)

completePD.to_excel('pubmed_treatment_info.xlsx', index=False)