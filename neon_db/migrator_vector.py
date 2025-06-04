import psycopg2
import ast

# Your Neon PostgreSQL connection string
conn = psycopg2.connect("postgresql://neondb_owner:npg_Jcn8LGTStZ3u@ep-still-hat-a66dlf3g-pooler.us-west-2.aws.neon.tech/neondb?sslmode=require")
cur = conn.cursor()

# Select text vectors and pmids from temp table
cur.execute("SELECT pmid, vector FROM temp.vectorexamples")
rows = cur.fetchall()

for pmid, vector_str in rows:
    try:
        # Safely convert text '[0.1, 0.2, ...]' to Python list
        vector = ast.literal_eval(vector_str)

        # Optional: validate length
        if len(vector) != 768:
            print(f"Skipping pmid {pmid}: invalid vector length")
            continue

        # Update pubmed_papers with real vector
        cur.execute("""
            UPDATE updated_treatment_data.pubmed_papers
            SET embedding = %s
            WHERE pmid = %s;
        """, (vector, pmid))

    except Exception as e:
        print(f"Error on pmid {pmid}: {e}")

conn.commit()
cur.close()
conn.close()
