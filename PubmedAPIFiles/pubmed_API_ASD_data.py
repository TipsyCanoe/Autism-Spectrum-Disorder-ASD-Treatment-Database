import pandas as pd
from Bio import Entrez
from pathlib import Path
from complete_query import importer
import os
from datetime import datetime

# Get all topics from df
def get_topics(topics):
    if topics:
        for i in range(len(topics)):
            if 'and' in topics[i]:
                # Remove 'and' for Pubmed
                topics[i] = topics[i].replace('and', '')
        topic_queries = ['{}'.format(topic) for topic in topics]
    return topic_queries

# Building full query
def build_ASD_query(author_queries, topic_queries):
    queries = []
    for i in range(max(len(author_queries), len(topic_queries))):
        author = author_queries
        if('Unknown' in author):
            author = ''
        topic = topic_queries[i]
        if('Unknown' in topic):
            topic = ''
        queries.append(str(author_queries[i]) + ' AND ' + str(topic) + ' AND ' + '(randomized controlled trial[Publication Type] OR "Clinical Trial"[Publication Type])')
    return queries

if __name__ == '__main__':
    # May need to adjust path depending on what directory you run this in
    path = str(Path.cwd()) + '/../PubmedAPIFiles/Harle - ASD bio tx .xlsx'

    # Null check
    if not os.path.isfile(path):
        print(f"File not found: {path}; skipping")
        exit(1)

    oldDf = pd.read_excel(path, 'Tx Data', usecols="A,D")

    # Getting file info
    oldDf.fillna("Unknown")
    authors = oldDf['First Author'].tolist() 
    topics = oldDf['Treatment name'].tolist()  

    topic_queries = get_topics(topics)
    
    queries = build_ASD_query(authors, topic_queries)

    completePD = pd.DataFrame()

    for query in queries:
        completePD = importer.importPapers(completePD, query)

    today = datetime.now().strftime("%Y-%m-%d")
    completePD.to_excel(f"pubmed_ASD_info_{today}.xlsx", index=False)
    

