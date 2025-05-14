CREATE SCHEMA IF NOT EXISTS updated_treatment_data;


CREATE TABLE updated_treatment_data.journals (

    id SERIAL PRIMARY KEY,

    journal VARCHAR(255) NOT NULL,

    publisher VARCHAR(255),

    impact_factor FLOAT,

    issn VARCHAR(20),

    url VARCHAR(255)

);


CREATE TABLE updated_treatment_data.treatments (

    id SERIAL PRIMARY KEY,

    "Treatment name" VARCHAR(255) UNIQUE NOT NULL,

    "Duration" VARCHAR(255),

    "Primary Outcome Area" VARCHAR(255),

    "Primary Outcome Measures" TEXT,

    "Results: Primary measure" TEXT,

    embedding VECTOR(768),

    journal_id INTEGER REFERENCES updated_treatment_data.journals(id)

);


CREATE TABLE updated_treatment_data.pubmed_papers (

    pmid VARCHAR(20) PRIMARY KEY,

    doi VARCHAR(255),

    title TEXT,

    pub_date DATE,

    abstract TEXT,

    authors TEXT,

    journal VARCHAR(255),

    keywords TEXT,

    url VARCHAR(255),

    affiliations TEXT,

    embedding VECTOR(768)

);


CREATE TABLE updated_treatment_data.treatment_pubmed_link (

    treatment_id INTEGER REFERENCES updated_treatment_data.treatments(id),

    pmid VARCHAR(20) REFERENCES updated_treatment_data.pubmed_papers(pmid),

    PRIMARY KEY (treatment_id, pmid)

);


DROP TABLE IF EXISTS pubmed_papers_info;
CREATE TABLE pubmed_papers_info (

    pmid VARCHAR(20),

    doi VARCHAR(255),

    title TEXT,

    pub_date DATE,

    abstract TEXT,

    authors TEXT,

    journal VARCHAR(255),

    keywords TEXT,

    url VARCHAR(255),

    affiliations TEXT,

    vector VECTOR(768)

);


DROP TABLE IF EXISTS vector_examples;
CREATE TABLE vector_examples (

    "Treatment name" VARCHAR(255),

    "Duration" VARCHAR(255),

    "Primary Outcome Area" VARCHAR(255),

    "Primary Outcome Measures" TEXT,

    "Results: Primary measure" TEXT,

    vector VECTOR(768)

);


INSERT INTO updated_treatment_data.journals (journal)

SELECT DISTINCT journal FROM pubmed_papers_info;


INSERT INTO updated_treatment_data.treatments (

    "Treatment name", "Duration", "Primary Outcome Area",

    "Primary Outcome Measures", "Results: Primary measure", embedding

)

SELECT DISTINCT

    "Treatment name", "Duration", "Primary Outcome Area",

    "Primary Outcome Measures", "Results: Primary measure", vector

FROM vector_examples;


INSERT INTO updated_treatment_data.pubmed_papers (

    pmid, doi, title, pub_date, abstract, authors,

    journal, keywords, url, affiliations, embedding

)

SELECT DISTINCT

    pmid, doi, title, pub_date, abstract, authors,

    journal, keywords, url, affiliations, vector

FROM pubmed_papers_info;


INSERT INTO updated_treatment_data.treatment_pubmed_link (treatment_id, pmid)

SELECT t.id, p.pmid

FROM updated_treatment_data.treatments t

JOIN updated_treatment_data.pubmed_papers p

ON t."Treatment name" = p.title;


SELECT * FROM updated_treatment_data.treatments LIMIT 5;

SELECT * FROM updated_treatment_data.pubmed_papers LIMIT 5;

SELECT * FROM updated_treatment_data.journals LIMIT 5;

SELECT * FROM updated_treatment_data.treatment_pubmed_link LIMIT 5;
