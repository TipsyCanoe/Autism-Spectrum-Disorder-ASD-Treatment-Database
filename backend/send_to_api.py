from openai import OpenAI

client = OpenAI()

batch_input_file = client.files.create(
    file=open("batch_requests.jsonl", "rb"),
    purpose="batch"
)

print(f"Uploaded file: {batch_input_file.filename}")
print(f"File ID: {batch_input_file.id}")

batch = client.batches.create(
    input_file_id=batch_input_file.id,
    endpoint="/v1/chat/completions",
    completion_window="24h",
    metadata={"description": "batch 1 test"}
)

print("Batch created successfully!")
print(f"Batch ID: {batch.id}")
print(f"Status: {batch.status}")
print(f"Created at: {batch.created_at}")