import pandas as pd
import json
from Bio import Entrez
from pathlib import Path

'''
- The main puller for the data
- Plan on implementing automatic pulling from a csv or text file 
depending on if the terms fit
- Also pipe to metapub maybe
- Ask me for the API key if you want to run it, or just make an account of your own and fill in your own info
- Github doesn't like it wihen keys are exposed
'''
# Set the email address to avoid any potential issues with Entrez
Entrez.email = 'loa4@wwu.edu'
#Entrez.api_key = ''

queries = []
topic_queries = []


# Will be adjusted
general_terms = ['Autism Spectrum Disorder[MeSH]', 'Autism[Title/Abstract]', 'Autistic Disorder[Title/Abstract]', 'ASD[Title/Abstract]']
biologics_terms = ['"Biological Products"[MeSH]', 'Biologic*[Title/Abstract]', '"Antibodies, monoclonal"[MeSH]', 'Cytokine Inhibitors[MeSH]', 'Immunotherapy[MeSH]']
treatment_terms = ['rTMS[Title/Abstract]', '"Transcranial Magnetic Stimulation"[MeSH]', 'ECT[Title/Abstract]', '"Electroconvulsive Therapy"[MeSH]', 
'Behavioral Therapy[MeSH]', '"Psychotherapy, Group"[MeSH]', 'Cognitive Therapy[MeSH]', 'Psychosocial Intervention[Title/Abstract]']
treatment_types_terms = ['randomized controlled trial[Publication Type]', 'randomized[Title/Abstract]', 'placebo[Title/Abstract]', '"Clinical Trial"[Publication Type]']

query_boxes = [general_terms, biologics_terms, treatment_terms, treatment_types_terms]

for terms in query_boxes:
    topic_queries = ['{}'.format(topic) for topic in terms]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

full_query = queries[0] + ' AND (' + queries[1] + ' OR ' + queries[2] + ') AND ' + queries[3]
print(full_query)

completePD = pd.DataFrame()

# Search PubMed for relevant records
handle = Entrez.esearch(db='pubmed', retmax=1000, term=full_query)
record = Entrez.read(handle)
id_list = record['IdList']

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['PMID', 'DOI', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

count = 0

# Fetch information for each record in the id_list
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(handle)

    # Process each PubMed article in the respons e
    for record in records['PubmedArticle']:
        if record.get('PubmedData'):
            if record['PubmedData'].get('ArticleIdList'):
                for other_id in record['PubmedData']['ArticleIdList']:
                    if 'doi' in other_id.attributes.values():
                        doi = (other_id.title()) 
        title = record['MedlineCitation']['Article']['ArticleTitle']
        abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in record['MedlineCitation']['Article'] and 'AbstractText' in record['MedlineCitation']['Article']['Abstract'] else ''
        authors = ', '.join(author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in record['MedlineCitation']['Article'].get('AuthorList', []))
        
        affiliations = []
        for author in record['MedlineCitation']['Article'].get('AuthorList', []):
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

        df = pd.concat([df, new_row], ignore_index=True).drop_duplicates()
        '''
        count += 1
        if(count % 100 == 0):
            print(count)
        '''
        

# Save DataFrame to an Excel file
df.to_excel('pubmed_papers_info.xlsx', index=False)

