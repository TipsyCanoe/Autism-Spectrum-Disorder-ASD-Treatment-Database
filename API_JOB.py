import subprocess
from pathlib import Path

def run_scripts(script_paths):
    """
    Runs a list of Python scripts.

    Args:
        script_paths: A list of strings, where each string is the path to a Python script.
    """

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

if __name__ == "__main__":
    scripts = []
    file_heads = ['pubmed_API_ASD_data.py', 'pubmed_API_data.py', 'pubmed_API_treatment']

    for file_head in file_heads:
        path = str(Path.cwd()) + '/PubmedAPIFiles/' + file_head
        scripts.append(path)

    run_scripts(scripts)