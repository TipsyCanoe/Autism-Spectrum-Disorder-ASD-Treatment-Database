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
Entrez.api_key = 'c301edeca095efe481ce5e2a727560444908'

# May need to adjust path depending on what directory you run this in
path = str(Path.cwd()) + '/Harle - ASD bio tx .xlsx'

oldDf = pd.read_excel(path, 'Tx Data', usecols="A,D")

authors = oldDf['First Author'].tolist() 
topics = oldDf['Treatment name'].tolist()  

queries = []
author_queries = []
topic_queries = []
monthDict = {
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
    "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
    "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
}

if authors:
    author_queries = ['{}'.format(author) for author in authors]

if topics:
    for treatment in topics:
        if 'and' in treatment:
            treatment.replace('and', '')
    topic_queries = ['{}'.format(topic) for topic in topics]

for i in range(len(author_queries)):
    queries.append(author_queries[i] + ' AND ' + topic_queries[i] + ' AND ' + '(randomized controlled trial[Publication Type] OR "Clinical Trial"[Publication Type])')

completePD = pd.DataFrame()

for query in queries:
    # Search PubMed for relevant records
    handle = Entrez.esearch(db='pubmed', retmax=10, term=query)
    record = Entrez.read(handle)
    id_list = record['IdList']

    # DataFrame to store the extracted data
    df = pd.DataFrame(columns=['PMID', 'DOI', 'Title', 'Publication Date', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

    # Fetch information for each record in the id_list
    for pmid in id_list:
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        records = Entrez.read(handle)

        # Process each PubMed article in the response
        for record in records['PubmedArticle']:
            doi = ''
            # Print the record in a formatted JSON style
            if record.get('PubmedData'):
                if record['PubmedData'].get('ArticleIdList'):
                    for other_id in record['PubmedData']['ArticleIdList']:
                        if 'doi' in other_id.attributes.values():
                            doi = (other_id.title())
            title = record['MedlineCitation']['Article']['ArticleTitle']
            pubDate = ''
            pubDate_section = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('PubDate', {})

            if 'Year' in pubDate_section:
                year = pubDate_section['Year']
                month = pubDate_section.get('Month', '01')
                if month != '01':
                    month = monthDict.get(month)
                day = pubDate_section.get('Day', '01')
                pubDate = f"{year}-{month}-{day}"
            elif 'MedlineDate' in pubDate_section:
                pubDate = pubDate_section['MedlineDate']  # Fallback to textual date
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
                'Publication Date': [pubDate],
                'Abstract': [abstract],
                'Authors': [authors],
                'Journal': [journal],
                'Keywords': [keywords],
                'URL': [url],
                'Affiliations': [affiliations]
            })

            df = pd.concat([df, new_row], ignore_index=True)
    completePD = pd.concat([completePD, df], ignore_index=True).drop_duplicates()
    # print("Finished " + query)

# Save DataFrame to an Excel file
completePD.to_excel('Adjusted_ASD_Sheet_Info_V2.xlsx', index=False)

