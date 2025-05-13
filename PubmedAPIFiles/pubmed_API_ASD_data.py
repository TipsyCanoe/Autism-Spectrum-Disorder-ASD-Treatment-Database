import pandas as pd
from Bio import Entrez
from pathlib import Path
from complete_query import importer

# Set the email address to avoid any potential issues with Entrez
Entrez.email = 'loa4@wwu.edu'
#Entrez.api_key = ''

# May need to adjust path depending on what directory you run this in
path = str(Path.cwd()) + '/Harle - ASD bio tx .xlsx'

oldDf = pd.read_excel(path, 'Tx Data', usecols="A,D")

# Getting file info
authors = oldDf['First Author'].tolist() 
topics = oldDf['Treatment name'].tolist()  

# Query initialization
queries = []
author_queries = []
topic_queries = []

# Get all authors...
if authors:
    author_queries = ['{}'.format(author) for author in authors]

# Get all topics...
if topics:
    for treatment in topics:
        if 'and' in treatment:
            treatment.replace('and', '')
    topic_queries = ['{}'.format(topic) for topic in topics]

# Building full query
for i in range(len(author_queries)):
    queries.append(author_queries[i] + ' AND ' + topic_queries[i] + ' AND ' + '(randomized controlled trial[Publication Type] OR "Clinical Trial"[Publication Type])')

completePD = pd.DataFrame()

for query in queries:
    importer.importPapers(completePD, query)

completePD.to_excel('Adjusted_ASD_Sheet_Info_V2.xlsx', index=False)

