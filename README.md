# Autism-Spectrum-Disorder-ASD-Treatment-Database

Vision: To enhance the mental health of individuals with Autism Spectrum Disorder (ASD)
and their families by synthesizing psychiatric treatment knowledge for healthcare
professionals, patients, and families

The goal of this project is to create a database and front-end representation for families and
healthcare professionals looking to know more about autism. Resources will be pulled from various
published studies and articles.

## Running All Project Servers Simultaneously

Two helper scripts are provided in the root directory to simplify starting and stopping all the necessary servers for this project:

* `start_all_servers.sh`: This script will attempt to start:
  * The frontend development server (typically on port 3000).
  * The backend Python query API (on port 5000).
  * The backend Node.js job API (on port 5001).
* `stop_all_servers.sh`: This script will attempt to find and stop the processes running on ports 3000, 5000, and 5001.

**How to Use:**

1. **Navigate to the Project Root:**
   Open your terminal and ensure you are in the root directory.
2. **Make the Scripts Executable (One-Time Setup):
   If you haven't done so already, you'll need to give the scripts execute permissions:

   ```bash
   chmod +x start_all_servers.sh
   chmod +x stop_all_servers.sh
   ```
3. **To Start All Servers:**
   Execute the start script:

   ```bash
   ./start_all_servers.sh
   ```

   The script will output the Process IDs (PIDs) of the servers it attempts to start. These servers will run in the background.
   This may render the terminal unusable the duration of the server uptime. To get back to being able to use commands,
   hit Ctrl+c
4. **To Stop All Servers:**
   Execute the stop script:

   ```bash
   ./stop_all_servers.sh
   ```

   This will attempt to terminate the processes listening on the specified ports.

**Important Considerations:**

* **Python Command:** The `start_all_servers.sh` script automatically detects whether to use `python3` or `python` based on what's available on your system.
* **`lsof` Dependency:** The `stop_all_servers.sh` script relies on the `lsof` command to find processes by port. This command is standard on most Linux systems. If it's missing, you might need to install it (e.g., `sudo apt install lsof` on Debian/Ubuntu-based systems).
* **Permissions:** If you encounter permission issues when stopping servers, you might need to run the `stop_all_servers.sh` script with `sudo` (e.g., `sudo ./stop_all_servers.sh`), especially if the servers were started with elevated privileges or by a different user.
* **Manual Closures:** If the frontend server (`npm start`) opens a new browser tab or terminal window, you may still need to close that manually after running `stop_all_servers.sh`.

If you want to run the frontend testing website locally, get into the frontend directory, and then the testing-website directory inside (``cd frontend/testing-website``)
Then, run ``npm install``. After that, running ``npm start`` should be sufficient to deploy locally. if you ctrl-z to
end the process in the terminal, the process still may be running on port 3000 (I think it's specified). run ``fuser -k 3000/tcp`` or ``npx kill-port 3000``
which should kill the process.

To run the dummy backend, temporarily, in the backend directory, run ``npm run dev``. In a different terminal, for the frontend, run the ``npm start`` command in the correct directory, as detailed above. The backend should run locally on port 5001.

For testing purposes, run ``npm test`` inside the frontend/testing-website/ directory. To visualize code coverage with the tests, run ``npm test -- --coverage --watchAll=false``

## PubMed API Extraction

Source: https://youtu.be/sGC66q45BX4

I started the code base from a public Github repository, referenced from the above Youtube video.
To increase the number of requests I can make per second to the API, I have included my email. I have not included my API key, as Github does not like exposed secrets. If you want to make faster requests, enter your own email and API key, freely obtainable through signing in with your account.

I have three extraction files currently; one based on previous data, one crafted based on common treatments and usages, and one based on lesser known treatments. The query based one is ``pubmed_API_data.py``, the common treatments is `pubmed_API_ASD_data.py`, and the lesser know treatments is `pubmed_API_treatment.py`. Running each will output an XML sheet of `pubmed_papers_info.xlsx`, `pubmed_ASD_info.xlsx`, and `pubmed_treatment_info.xlsx` respectively.

`convert_excel_to_json.py` is as the name implies; it converts the excel spreadsheet to a JSON. Currently, it expects the three excel sheets listed above.

Run `API_JOB.py` to start running the job automatically manually. If the website is started, it will update the database if "Update Database" is pressed, or every Monday at 2 am. Change the time in /backend/scheduler.js, line 65. Use crontab.guru to help figure out an exact time if unfamiliar with cron. 

Dependencies (if using Python 3, just add 3 after pip or python): ``sh pip install pandas biopython``

### Using the medBERT/LLM

1. Create a virtual environment if it's your first time running
   python3 -m venv venv
2. Activate the environment
   source venv/bin/activate
3. Install the packages
   pip install -r requirements.txt
4. You should now be able to run the LLM and MedBERT

### Using the backend

1. Create a virtual environment if it's your first time running
   python3 -m venv venv
2. Activate the environment
   source venv/bin/activate
3. Install the packages (in backend if you don't want it to take a while)
   pip install -r requirements.txt
4. run app.py and the frontend on seperate terminals
