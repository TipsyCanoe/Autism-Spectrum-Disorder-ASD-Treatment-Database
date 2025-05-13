import pandas as pd
import json
import re
from Bio import Entrez
from pathlib import Path

monthDict = {
        "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
        "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
        "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
    }

class importer:
    Entrez.email = 'loa4@wwu.edu'
    Entrez.api_key = 'c301edeca095efe481ce5e2a727560444908'

    def importPapers(completePD, query):
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
        completePD = pd.concat([completePD, df], ignore_index=True).drop_duplicates()
        return completePD