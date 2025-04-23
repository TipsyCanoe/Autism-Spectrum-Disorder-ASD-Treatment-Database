import pandas as pd
import json
import re
from Bio import Entrez
from pathlib import Path

def expand(text):
    pattern = r"([^()\sANDOR]+(?: [^()\sANDOR]+)*)(?=\s*(?:AND|OR|\(|\)|$))|([()])|(AND|OR)"
    #tokens = re.findall(r'\w+|[()]', text)    
    tokens = re.findall(pattern, text)
    # print(tokens)
    result = []

    for token in tokens:
        if token in {'AND', 'OR', '(', ')'}:
            result.append(token)
        else:
            term = token.strip()
            # print(term)
            if term:
                result.append(f'("{term}"[MeSH terms] OR "{term}"[All Fields])')

    return ' '.join(result)

monthDict = {
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
    "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
    "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
}

'''
- The main puller for the data
- Github doesn't like it wihen keys are exposed
'''
# Set the email address to avoid any potential issues with Entrez
Entrez.email = 'loa4@wwu.edu'
#Entrez.api_key = ''

# May need to adjust path depending on what directory you run this in
path = str(Path.cwd()) + '/Treatment_Names.xlsx'

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
            # treatment_new = expand(treatments[i]);
            # print(treatment_new)
    treatment_queries = ['({})'.format(treatment) for treatment in treatments]
# print(treatment_queries)


general_terms = ['Autism Spectrum Disorder[MeSH]', 'Autism[Title/Abstract]', 'Autistic Disorder[Title/Abstract]', 'ASD[Title/Abstract]']

treatment_types_terms = ['randomized controlled trial[Publication Type]', 'randomized[Title/Abstract]', 'placebo[Title/Abstract]', '"Clinical Trial"[Publication Type]']

for terms in [general_terms, treatment_types_terms]:
    topic_queries = ['{}'.format(topic) for topic in terms]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

for treatment in treatment_queries:
    full_queries.append(treatment + ' AND ' + queries[0] + ' AND ' + queries[1])
#print(full_queries[0])

completePD = pd.DataFrame()

count = 0

for query in full_queries:
    # Search PubMed for relevant records
    handle = Entrez.esearch(db='pubmed', retmax=1000, term=query)
    record = Entrez.read(handle)
    id_list = record['IdList']

    # DataFrame to store the extracted data
    df = pd.DataFrame(columns=['pmid', 'doi', 'title', 'pub_date', 'abstract', 'authors', 'journal', 'keywords', 'url', 'affiliations'])

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
            
            # Safely handle 'AuthorList'
            authors = ''
            affiliations = []

            if 'AuthorList' in record['MedlineCitation']['Article']:
                authors = ', '.join(
                    author.get('LastName', '') + ' ' + author.get('ForeName', '')
                    for author in record['MedlineCitation']['Article']['AuthorList']
                )
                for author in record['MedlineCitation']['Article']['AuthorList']:
                    if 'AffiliationInfo' in author and author['AffiliationInfo']:
                        affiliations.append(author['AffiliationInfo'][0]['Affiliation'])

            affiliations = '; '.join(set(affiliations))

            journal = record['MedlineCitation']['Article']['Journal']['Title']
            keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''
            url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

            new_row = pd.DataFrame({
                'pmid': [pmid],
                'doi': [doi],
                'title': [title],
                'pub_date': [pubDate],
                'abstract': [abstract],
                'authors': [authors],
                'journal': [journal],
                'keywords': [keywords],
                'url': [url],
                'affiliations': [affiliations]
            })

            df = pd.concat([df, new_row], ignore_index=True)
            count += 1
            if(count % 100 == 0):
                print(count)
    completePD = pd.concat([completePD, df], ignore_index=True).drop_duplicates()
    

# Save DataFrame to an Excel file
completePD.to_excel('Individual_Treatment_Queries.xlsx', index=False)

