import requests
from metapub import PubMedFetcher
from metapub.findit import FindIt
from pypdf import PdfReader

# Dependencies: pypdf and metapub
# Run program and it should pull the paper and download the pdf and print the text
# of the first page into a txt file
def main():
    fetch = PubMedFetcher()
    # article = fetch.article_by_pmid('29732299')
    src = FindIt(doi='10.1038/npp.2016.237')
    url = src.url
    response = requests.get(url)
    file_Path = 'test_paper.pdf'

    if response.status_code == 200:
        with open(file_Path, 'wb') as file:
            file.write(response.content)
        print('File downloaded successfully')
    else:
        print('Failed to download file')

    reader = PdfReader("test_paper.pdf")
    page = reader.pages[0]
    text = page.extract_text()
    file = open("test_paper_pg_1.txt", "w")
    file.write(text)
    file.close()
    

if __name__ == "__main__":
    main()