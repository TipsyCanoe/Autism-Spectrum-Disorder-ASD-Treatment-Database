# Updating the Neon Database from WSL

This guide explains how to connect to your [Neon PostgreSQL database](https://neon.tech/) from **WSL (Windows Subsystem for Linux)** and load updated data using `.csv` files.

---

## Prerequisites

1. You have [psql](https://www.postgresql.org/docs/current/app-psql.html) installed in your WSL environment.
2. Your Neon database connection string (Postgres URI) is ready (from the Neon Console).
3. All `.csv` files are cleaned and correctly formatted.
4. You're using **PostgreSQL v14+** to match Neon server requirements.

---

## Step 1: Set Your Connection String

You can connect using either:

### Option A: Direct psql command

```bash
psql "postgres://<user>:<password>@<host>/<dbname>?sslmode=require"
```

Replace:

* `<user>`, `<password>`, `<host>`, and `<dbname>` with your values from Neon.

Or, copy and past the connection string from the neondb dashboard:

```bash
psql "connection_string"
```

---

## Step 2: Ensure File Paths Are Linux-Readable

If your `.csv` files are in a Windows directory, copy them into your Linux home directory for compatibility:

```bash
cp /mnt/c/Users/yourname/Path/to/file.csv ~/file.csv
```

---

## Step 3: Import Data Using `\copy`

Once inside `psql`, use the `\copy` command (not `COPY`) so PostgreSQL reads the file **from the client (your WSL machine)**:

### Example – Importing Journals:

```sql
\copy updated_treatment_data.journals(journal, publisher, impact_factor, issn, url)
FROM '/home/youruser/cleaned_journals.csv' DELIMITER ',' CSV HEADER;
```

### Example – Importing PubMed Papers:

Make sure your `pub_date` is properly formatted as `YYYY-MM-DD`. Then:

```sql
\copy updated_treatment_data.pubmed_papers(pmid, doi, title, pub_date, abstract, authors, journal, keywords, url, affiliations)
FROM '/home/youruser/cleaned_pubmed_papers.csv' DELIMITER ',' CSV HEADER;
```

---

## 🧪 Troubleshooting

* **SSL SYSCALL error / EOF detected**: This may occur if the file path is invalid or the connection drops — copy files to WSL (`~/`) and retry.
* **Date Format Error**: Ensure your dates are in `YYYY-MM-DD` format (consider using a script like `pub_date_clean.py` to fix malformed dates).
* **Value too long for type**: Check column definitions (`VARCHAR(n)`) and adjust schema using `ALTER TABLE`.

---

## Optional: Clean Up

If you've used temporary copies, remove them:

```bash
rm ~/cleaned_journals.csv ~/cleaned_pubmed_papers.csv
```

---

