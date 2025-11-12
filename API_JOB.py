import subprocess
from pathlib import Path
from datetime import datetime
import os
from dotenv import load_dotenv

# runs the designated series of scripts
def run_scripts(script_paths):
    for script_path in script_paths:
        try:
            print("Running " + script_path)
            result = subprocess.run(['python3', script_path], capture_output=True, text=True, check=True)
            if result.stderr:
                print(f"--- Errors from {script_path} ---")
                print(result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"--- Error running {script_path} ---")
            print(e.stderr)
        except FileNotFoundError:
            print(f"Error: Python interpreter not found.")
        except Exception as e:
            print(f"An unexpected error occurred while running {script_path}: {e}")
    
def update_pull_data():
    today = datetime.now().strftime("%Y-%m-%d")
    time = datetime.now().strftime("%H:%M:%S")
    with open('.env', 'w') as f:
        f.write(f"LAST_PULL_DATE={today}\n")
        f.write(f"LAST_PULL_TIME={time}\n")

if __name__ == "__main__":
    scripts = []
    file_heads = ['pubmed_API_ASD_data.py', 'pubmed_API_data.py', 'pubmed_API_treatment.py', 'convert_excels_to_csvs.py']

    for file_head in file_heads:
        path = str(Path.cwd()) + '/PubmedAPIFiles/' + file_head
        scripts.append(path)

    load_dotenv()
    run_scripts(scripts)
    update_pull_data()