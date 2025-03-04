import pandas as pd
import json
from Bio import Entrez
from pathlib import Path
'''
Things I removed and should maybe remember:
- Full queries
- [Title/Abstract] and [Author]
'''
# Set the email address to avoid any potential issues with Entrez
Entrez.email = 'loa4@wwu.edu'
#Entrez.api_key = ''

# May need to adjust path depending on what directory you run this in
path = str(Path.cwd()) + '/PubmedAPIFiles/Harle - ASD bio tx .xlsx'
#print(path)

oldDf = pd.read_excel(path, 'Tx Data', usecols="A,D")

# Define lists of authors and topics

authors = oldDf['First Author'].tolist() # Example authors, adjust as needed
topics = oldDf['Treatment name'].tolist()  # Example topics, adjust as needed

# Define date range
#date_range = '("2012/03/01"[Date - Create] : "2022/12/31"[Date - Create])'

# Build the query dynamically based on the available authors and topics
queries = []
author_queries = []
topic_queries = []

if authors:
    #author_queries = ['{}[Author]'.format(author) for author in authors]
    author_queries = ['{}'.format(author) for author in authors]
    #print(author_queries)
    #queries.append('(' + ' OR '.join(author_queries) + ')')

if topics:
    for treatment in topics:
        if 'and' in treatment:
            treatment.replace('and', '')
    topic_queries = ['{}'.format(topic) for topic in topics]
    #queries.append('(' + ' OR '.join(topic_queries) + ')')

for i in range(len(author_queries)):
    queries.append(author_queries[i] + ' AND ' + topic_queries[i])

'''
print(author_queries)
print(topic_queries)
print(queries)
'''
completePD = pd.DataFrame()

for query in queries:
    # Search PubMed for relevant records
    handle = Entrez.esearch(db='pubmed', retmax=10, term=query)
    record = Entrez.read(handle)
    id_list = record['IdList']

    # DataFrame to store the extracted data
    df = pd.DataFrame(columns=['PMID', 'DOI', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

    # Fetch information for each record in the id_list
    for pmid in id_list:
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        records = Entrez.read(handle)

        # Process each PubMed article in the response
        for record in records['PubmedArticle']:
            # Print the record in a formatted JSON style
            #print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime
            if record.get('PubmedData'):
                if record['PubmedData'].get('ArticleIdList'):
                    for other_id in record['PubmedData']['ArticleIdList']:
                        if 'doi' in other_id.attributes.values():
                            doi = (other_id.title())
            title = record['MedlineCitation']['Article']['ArticleTitle']
            abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in record['MedlineCitation']['Article'] and 'AbstractText' in record['MedlineCitation']['Article']['Abstract'] else ''
            authors = ', '.join(author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in record['MedlineCitation']['Article']['AuthorList'])
            
            affiliations = []
            for author in record['MedlineCitation']['Article']['AuthorList']:
                if 'AffiliationInfo' in author and author['AffiliationInfo']:
                    affiliations.append(author['AffiliationInfo'][0]['Affiliation'])
            affiliations = '; '.join(set(affiliations))

            journal = record['MedlineCitation']['Article']['Journal']['Title']
            keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''
            url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

            new_row = pd.DataFrame({
                'PMID': [pmid],
                'DOI': [doi],
                'Title': [title],
                'Abstract': [abstract],
                'Authors': [authors],
                'Journal': [journal],
                'Keywords': [keywords],
                'URL': [url],
                'Affiliations': [affiliations]
            })

            df = pd.concat([df, new_row], ignore_index=True)
    completePD = pd.concat([completePD, df], ignore_index=True).drop_duplicates()
    print("Finished " + query)

# Save DataFrame to an Excel file
completePD.to_excel('Adjusted_ASD_Sheet_Info.xlsx', index=False)

