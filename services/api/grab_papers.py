import pandas as pd
from sqlalchemy import create_engine

engine = create_engine('postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require&channel_binding=require')

df = pd.read_sql("SELECT pmid, doi, title, pub_date, abstract, authors, journal, keywords, url, affiliations, embedding FROM updated_treatment_data.pubmed_papers", engine)

df.to_csv('grabbed_papers.csv', index=False)