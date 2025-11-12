from pubmed_API_ASD_data import get_topics, build_ASD_query
from pubmed_API_data import get_full_query, get_terms
from pubmed_API_treatment import get_full_treatment_queries, get_treatments_queries
from complete_query import importer
import shutil
from pathlib import Path
import pandas as pd

def test_get_topics():
    topics = ['ariprizole', 'health and 9']
    topic_query = get_topics(topics)
    assert topic_query == ['ariprizole', 'health  9']

    topics = ['tdms', 'autism spectrum disorder']
    topic_query = get_topics(topics)
    assert topic_query == topics

def test_build_ASD_query():
    authors = ['John Doe', 'Unknown']
    topics = ['ariprizole', 'health and 9']
    topic_query = get_topics(topics)

    asd_query = build_ASD_query(authors, topic_query)
    assert asd_query == ['John Doe AND ariprizole AND (randomized controlled trial[Publication Type] OR "Clinical Trial"[Publication Type])', 
                         'Unknown AND health  9 AND (randomized controlled trial[Publication Type] OR "Clinical Trial"[Publication Type])']
    
def test_get_full_query(tmp_path, monkeypatch):
    file_cp = str(Path.cwd()) + '/biologic_terms.txt'
    file_cp2 = str(Path.cwd()) + '/social_terms.txt'
    monkeypatch.chdir(tmp_path)
    temp_file = tmp_path / "biologic_terms.txt"
    temp_file2 = tmp_path / "social_terms.txt"
    shutil.copy(file_cp, temp_file)
    shutil.copy(file_cp2, temp_file2)
    query = '(Autism Spectrum Disorder[MeSH] OR Autism[Title/Abstract] OR Autistic Disorder[Title/Abstract] OR ASD[Title/Abstract]) AND (("Biological Products"[MeSH] OR Biologic*[Title/Abstract] OR "Antibodies, monoclonal"[MeSH] OR \'Cytokine Inhibitors[MeSH]\' OR \'Immunotherapy[MeSH]\') OR (rTMS[Title/Abstract] OR "Transcranial Magnetic Stimulation"[MeSH] OR ECT[Title/Abstract] OR "Electroconvulsive Therapy"[MeSH] OR Behavioral Therapy[MeSH] OR "Psychotherapy, Group"[MeSH] OR Cognitive Therapy[MeSH] OR Psychosocial Intervention[Title/Abstract])) AND (randomized controlled trial[Publication Type] OR randomized[Title/Abstract] OR placebo[Title/Abstract] OR "Clinical Trial"[Publication Type])'

    # Checking against known files
    assert get_full_query() == query

def test_get_terms(tmp_path, monkeypatch):
    default = ['testing...']

    file_cp = str(Path.cwd()) + '/biologic_terms.txt'
    monkeypatch.chdir(tmp_path)
    temp_file = tmp_path / "biologic_terms.txt"
    shutil.copy(file_cp, temp_file)
    assert get_terms("biologic_terms.txt", default) == ['"Biological Products"[MeSH]', 'Biologic*[Title/Abstract]', '"Antibodies, monoclonal"[MeSH]',"'Cytokine Inhibitors[MeSH]'", "'Immunotherapy[MeSH]'"]
    
    # No existence check
    assert get_terms("hello.txt", default) == default

def test_full_treatment_queries():
    treatments = ['normal', 'test And', 'test Or Clause', 'test aND and oR clause']
    treatment_query = get_treatments_queries(treatments)
    queries = ['test1', 'test2']

    assert get_full_treatment_queries(treatment_query, queries) == ['(normal) AND test1 AND test2', '(test AND) AND test1 AND test2', '(test OR Clause) AND test1 AND test2', '(test AND AND OR clause) AND test1 AND test2']

def test_get_treatment_queries():
    treatments = ['normal', 'test And', 'test Or Clause', 'test aND and oR clause']
    treatment_query = get_treatments_queries(treatments)
    assert treatment_query == ['(normal)', '(test AND)', '(test OR Clause)', '(test AND AND OR clause)']

def test_importer_part():
    df = pd.DataFrame()

    columns = ['pmid', 'doi', 'title', 'pub_date', 'abstract', 'authors', 'journal', 'keywords', 'url', 'affiliations']

    df = importer.importPapers(df, "ajsdklfjajsdfajsdiofj")
    assert list(df.columns) == columns
