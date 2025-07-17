# Autism-Spectrum-Disorder-ASD-Treatment-Database

Vision: To enhance the mental health of individuals with Autism Spectrum Disorder (ASD) and their families by synthesizing psychiatric treatment knowledge for healthcare professionals, patients, and families.

## Prerequisites

Before you begin, make sure you have the following installed on your system:

- **Python 3.8+**
  - [Download for Windows/macOS/Linux](https://www.python.org/downloads/)
  - On Linux:
    ```bash
    sudo apt update
    sudo apt install python3 python3-venv python3-pip
    ```
  - On macOS (with Homebrew):
    ```bash
    brew install python
    ```
  - On Windows: Download and run the installer from the Python website above.

- **Node.js & npm**
  - [Download for all platforms](https://nodejs.org/)

---

## Quick Start Guide

### 1. First-Time Setup

**A. Python Backend**
1. Open a terminal in the project root after cloning the project locally.
2. Create and activate a Python virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. Install backend dependencies:
   ```bash
   pip install -r backend/requirements.txt
   ```

**B. Frontend**
1. Go to the frontend directory:
   ```bash
   cd frontend/testing-website
   ```
2. Install frontend dependencies:
   ```bash
   npm install
   ```

### 2. Running the Project

**A. Start All Servers (from project root):**
```bash
chmod +x start_all_servers.sh stop_all_servers.sh  # One-time
./start_all_servers.sh
```
- This starts:
  - Frontend (port 3000)
  - Python backend (port 5000)
  - Node.js job backend (port 5001)

**B. Stop All Servers:**
```bash
./stop_all_servers.sh
```

### 3. Manual Server Control

- **Frontend only:**  
  ```bash
  cd frontend/testing-website
  npm start
  ```
- **Python backend only:**  
  ```bash
  cd backend
  source ../venv/bin/activate
  python app.py
  ```
- **Node.js backend only:**  
  ```bash
  cd backend
  npm run dev
  ```

### 4. Testing

- **All tests:**  
  ```bash
  ./run_all_tests.sh
  ```
- **Frontend tests:**  
  ```bash
  cd frontend/testing-website
  npm test
  ```
- **Backend tests:**  
  ```bash
  cd backend/tests
  ./run_tests.sh
  ```

---

## PubMed API Extraction

- See `pubmed_API_data.py`, `pubmed_API_ASD_data.py`, and `pubmed_API_treatment.py` for extraction scripts.
- Convert Excel to JSON with `convert_excel_to_json.py`.
- Run `API_JOB.py` to manually update the database or use the website's update feature.
- Change scheduled update time in `/backend/scheduler.js` (see [crontab.guru](https://crontab.guru)).

## Using medBERT/LLM

1. Create a virtual environment (see above).
2. Activate the environment.
3. Install packages:
   ```bash
   pip install -r requirements.txt
   ```
4. You can now run the LLM and MedBERT scripts.

---

For more details, see comments in the code and scripts in the repository.
