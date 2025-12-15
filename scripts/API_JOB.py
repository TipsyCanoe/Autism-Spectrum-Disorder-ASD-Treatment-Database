import subprocess
from pathlib import Path
from datetime import datetime
import os
import sys
from dotenv import load_dotenv
import urllib.request
import urllib.error

# runs the designated series of scripts
def run_scripts(script_paths):
    for script_path in script_paths:
        try:
            print("Running " + script_path)
            result = subprocess.run([sys.executable, script_path], capture_output=True, text=True, check=True)
            if result.stderr:
                print(f"--- Errors from {script_path} ---")
                print(result.stderr)
            else:
                print(f"--- Finished {script_path} ---")
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
    
    env_path = Path.cwd() / '.env'
    env_lines = []
    if env_path.exists():
        with open(env_path, 'r') as f:
            env_lines = f.readlines()
    
    new_lines = []
    keys_updated = {'LAST_PULL_DATE': False, 'LAST_PULL_TIME': False}
    
    for line in env_lines:
        if line.startswith('LAST_PULL_DATE='):
            new_lines.append(f"LAST_PULL_DATE={today}\n")
            keys_updated['LAST_PULL_DATE'] = True
        elif line.startswith('LAST_PULL_TIME='):
            new_lines.append(f"LAST_PULL_TIME={time}\n")
            keys_updated['LAST_PULL_TIME'] = True
        else:
            new_lines.append(line)
            
    if not keys_updated['LAST_PULL_DATE']:
        new_lines.append(f"LAST_PULL_DATE={today}\n")
    if not keys_updated['LAST_PULL_TIME']:
        new_lines.append(f"LAST_PULL_TIME={time}\n")
        
    with open(env_path, 'w') as f:
        f.writelines(new_lines)

def clear_cache():
    port = os.getenv('PYTHON_BACKEND_PORT', 5000)
    url = f"http://localhost:{port}/api/clear-cache"
    try:
        req = urllib.request.Request(url, method='POST')
        with urllib.request.urlopen(req) as response:
            if response.status == 200:
                print("Successfully cleared backend cache.")
            else:
                print(f"Failed to clear cache. Status code: {response.status}")
    except urllib.error.URLError as e:
        print(f"Could not connect to backend to clear cache: {e}")
    except Exception as e:
        print(f"Error clearing cache: {e}")

if __name__ == "__main__":
    # Load variables from environment
    load_dotenv(override=True)
    today = datetime.now().strftime("%Y-%m-%d")
    last_pull_date = os.getenv('LAST_PULL_DATE')
    
    # Only pull if we haven't pulled yet today
    if last_pull_date != today:
        # Appending file paths for API scraping scripts
        scripts = []
        file_heads = [
            "pubmed_API_ASD_data.py",
            "pubmed_API_data.py",
            "pubmed_API_treatment.py",
            "convert_excels_to_csvs.py"
        ]
        for file_head in file_heads:
            path = str(Path.cwd()) + '/pubmed_api/' + file_head
            scripts.append(path)

        run_scripts(scripts)
        update_pull_data()
        print("Finished running all API scripts")

        # Running upload pipeline scripts
        path_to_uploaders = [str(Path.cwd()) + '/services/api/LLMPipeline.py', 
                       str(Path.cwd()) + '/neon_db/automated_csv_uploader.py']
        run_scripts(path_to_uploaders)
        
        # Clear the cache after all updates are done
        clear_cache()
    else:
        print("Already pulled today. Please try again tommorrow.")
    