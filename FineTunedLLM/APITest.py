from openai import OpenAI

# Initialize client
client = OpenAI()

# Upload your batch input file
batch_input_file = client.files.create(
    file=open("batchinput.jsonl", "rb"),
    purpose="batch"
)

print(f"Uploaded file: {batch_input_file.filename}")
print(f"File ID: {batch_input_file.id}")

# Create the batch job
batch = client.batches.create(
    input_file_id=batch_input_file.id,
    endpoint="/v1/chat/completions",
    completion_window="24h",
    metadata={"description": "batch 1 test"}
)

# Print info for confirmation
print("Batch created successfully!")
print(f"Batch ID: {batch.id}")
print(f"Status: {batch.status}")
print(f"Created at: {batch.created_at}")