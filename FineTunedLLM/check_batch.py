import sys
import os
from openai import OpenAI

def main():
    if len(sys.argv) != 2:
        print("Usage: python check_batch.py <batch_id>")
        sys.exit(1)

    batch_id = sys.argv[1]
    client = OpenAI()

    try:
        # Retrieve the batch info
        batch = client.batches.retrieve(batch_id)
        print(f"Batch ID: {batch.id}")
        print(f"Status: {batch.status}")
        print(f"Created at: {batch.created_at}")

        # If it's completed, download the output file
        if batch.status == "completed":
            if not hasattr(batch, "output_file_id") or batch.output_file_id is None:
                print("No output file found for this batch.")
                return

            output_file_id = batch.output_file_id
            print(f"Downloading output file: {output_file_id}")

            # Retrieve file content
            output_file = client.files.content(output_file_id)
            output_data = output_file.text

            # Ensure output directory exists
            output_dir = "outputs"
            os.makedirs(output_dir, exist_ok=True)

            # Save to local file
            output_filename = os.path.join(output_dir, f"{batch_id}_output.jsonl")
            with open(output_filename, "w", encoding="utf-8") as f:
                f.write(output_data)

            print(f"Output saved to {output_filename}")

        else:
            print("Batch not completed yet, try again later.")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()