# Neon Database CSV Upload Tool

Interactive Python tool for uploading CSV files to the Neon PostgreSQL database with validation, preview, and progress tracking.

---

## Features

- ✅ Interactive command-line interface with smart defaults
- ✅ CSV preview before upload (first 5 rows)
- ✅ File validation and error checking
- ✅ Automatic column name cleaning (lowercase, underscores)
- ✅ NULL value handling
- ✅ Upload progress tracking with duration
- ✅ Pre-configured for `jim_data.data_embedded` table

---

## Quick Start

### Installation

```bash
# Install required packages
pip3 install pandas sqlalchemy psycopg2-binary
```

Or use the requirements file:

```bash
pip3 install -r neon_db/requirements-interactive.txt
```

### Basic Usage

```bash
# Ensure Correct Directory
cd /path/to/Autism-Spectrum-Disorder-ASD-Treatment-Database

# Run the upload script
python3 neon_db/interactive_csv_upload.py
```

---

## How It Works

### Step 1: Run the Script

```bash
python3 neon_db/interactive_csv_upload.py
```

### Step 2: Enter CSV Filename

When prompted, enter the filename from `data/db_imports/` directory:

```
Enter CSV filename (in data/db_imports/): your_file.csv
```

**Note:** The script automatically prepends `data/db_imports/` to your filename.

**Examples:**
- `test.csv` - Test file with sample data
- `pubmed_papers_info.csv` - PubMed paper data
- `study_details_final.csv` - Study details data

### Step 3: Configure Upload (Optional)

Press Enter to accept defaults or customize:

```
Enter schema name [jim_data]: 
Enter table name [data_embedded]: 
Upload mode (append/replace) [append]: 
```

### Step 4: Review and Confirm

Review the CSV preview and upload summary, then confirm:

```
Proceed with upload? (y/N): y
```

---

## Example Session

```
🚀 Neon Database CSV Uploader
========================================
Enter CSV filename (in PubmedAPIFiles/): your_data.csv
Enter schema name [jim_data]: 
Enter table name [data_embedded]: 
Upload mode (append/replace) [append]: 

📊 CSV Preview (5 rows shown):
     pmid              doi                          title  pub_date
  1234567  10.1234/example  Study on ASD Treatment X...       2024-01-15
  2345678  10.2345/example  Behavioral Intervention...        2024-02-20
  ...

Columns: ['pmid', 'doi', 'title', 'pub_date', 'abstract', ...]
File size: 2.45 MB

🎯 Upload Summary:
   File: PubmedAPIFiles/your_data.csv
   Target: jim_data.data_embedded
   Mode: append

Proceed with upload? (y/N): y

🔄 Starting upload...
   Loaded 1234 rows, 30 columns
✅ Upload successful!
   Rows uploaded: 1234
   Duration: 3.21 seconds
   Table: jim_data.data_embedded
```

---

## Configuration

### Default Settings

| Setting | Default Value | Description |
|---------|--------------|-------------|
| Schema | `jim_data` | Database schema |
| Table | `data_embedded` | Target table |
| Mode | `append` | Upload mode |
| File Location | `data/db_imports/` | CSV directory |

### Upload Modes

**Append Mode (Default)**
- Adds new rows to existing table
- Preserves all existing data
- Safe for incremental updates

**Replace Mode**
- ⚠️ Deletes all existing table data
- Uploads fresh data
- Use with caution!

---

## CSV File Requirements

### Location
Place CSV files in the `data/db_imports/` directory at project root.

**Example CSV Locations:**
```
Autism-Spectrum-Disorder-ASD-Treatment-Database/
└── data/
    └── db_imports/
        ├── test.csv                    (example test file)
        ├── pubmed_papers_info.csv      (PubMed data)
        ├── study_details_final.csv     (study details)
        └── your_data.csv               (any CSV file)
```

### Format Requirements

**Column Names:** Automatically cleaned to lowercase with underscores
- `"PMID"` → `"pmid"`
- `"Pub Date"` → `"pub_date"`

**Vector Columns:** Must contain exactly **768 dimensions**
```csv
"embedding"
"[0.001,0.002,0.003,...,0.768]"
```

**Date Format:** `YYYY-MM-DD`
```csv
"pub_date"
"2024-01-15"
```

**NULL Values:** Empty cells automatically handled as SQL `NULL`

---

## Troubleshooting

### Common Errors

**File Not Found**
```
❌ File not found: data/db_imports/yourfile.csv
```
- Verify file exists in `data/db_imports/` directory
- Check spelling and `.csv` extension
- List files: `ls data/db_imports/`

**Vector Dimension Error**
```
❌ Upload failed: expected 768 dimensions, not 3
```
- Ensure `embedding` and `vector` columns have 768 values
- Use test data generator for proper format

**Duplicate Key Error**
```
❌ Upload failed: duplicate key value violates unique constraint
```
- Use `replace` mode to overwrite
- Or delete conflicting rows first

**Connection Error**
```
❌ Upload failed: could not connect to server
```
- Check internet connection
- Verify Neon database is active

---

## Database Operations

### Query Data

```sql
-- Search by PMID (example with test data)
SELECT pmid, title, treatment_name 
FROM jim_data.data_embedded
WHERE pmid = 99999999;

-- Search for real data by PMID
SELECT pmid, title, treatment_name 
FROM jim_data.data_embedded
WHERE pmid = 12345678;

-- Count all records
SELECT COUNT(*) FROM jim_data.data_embedded;

-- View recent uploads
SELECT pmid, title, pub_date
FROM jim_data.data_embedded
ORDER BY pub_date DESC
LIMIT 10;
```

### Delete Records

```sql
-- Delete specific record by PMID
DELETE FROM jim_data.data_embedded WHERE pmid = 99999999;

-- Delete all test records
DELETE FROM jim_data.data_embedded WHERE title LIKE 'TEST:%';

-- Delete multiple specific records
DELETE FROM jim_data.data_embedded 
WHERE pmid IN (99999999, 88888888, 77777777);
```

---

## Testing

### Using the Provided Test File

A sample test CSV with proper 768-dimensional vectors is provided for testing purposes:
```
data/db_imports/test.csv
```

This test file contains dummy data with PMID 99999999 to help you verify the upload process works correctly.

**Upload test data:**
```bash
python3 neon_db/interactive_csv_upload.py
# Enter: test.csv
# Press Enter for defaults
# Confirm with 'y'
```

**Verify upload:**
```sql
SELECT * FROM jim_data.data_embedded WHERE pmid = 99999999;
```

**Clean up:**
```sql
DELETE FROM jim_data.data_embedded WHERE pmid = 99999999;
```

### Uploading Your Own Data

To upload your actual data files:

1. Place your CSV in `data/db_imports/` directory
2. Ensure it matches the required format (see CSV File Requirements)
3. Run the script and enter your filename
4. Review the preview and confirm

**Example:**
```bash
python3 neon_db/interactive_csv_upload.py
# Enter: pubmed_papers_info.csv
# Press Enter for defaults
# Review preview and confirm with 'y'
```

---

## Advanced

### Custom Database Connection

Set environment variable:
```bash
export DATABASE_URL="postgresql://user:password@host/db?sslmode=require"
python3 neon_db/interactive_csv_upload.py
```

### For Large Files (>1GB)

Use the automated uploader with chunked processing:
```bash
python3 neon_db/automated_csv_uploader.py
```

### Manual psql Upload

For advanced users:
```bash
psql "connection_string"
\copy jim_data.data_embedded FROM 'data/db_imports/file.csv' DELIMITER ',' CSV HEADER;
```

---

## Resources

- [Neon Documentation](https://neon.tech/docs)
- [Pandas CSV Docs](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html)
- [SQLAlchemy Docs](https://docs.sqlalchemy.org/)

---

