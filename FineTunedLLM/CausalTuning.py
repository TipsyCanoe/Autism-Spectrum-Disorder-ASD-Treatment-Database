import torch
from transformers import AutoModelForCausalLM, AutoTokenizer, TrainingArguments, Trainer, DataCollatorForLanguageModeling
from peft import LoraConfig, get_peft_model, TaskType, prepare_model_for_kbit_training
from datasets import load_dataset
import os

model_name = "meta-llama/Llama-3.2-1B"
output_dir = "./fine_tuned_model"
data_file = "V3PR.json"


dataset = load_dataset("json", data_files=data_file)
dataset = dataset["train"].train_test_split(test_size=0.1)


tokenizer = AutoTokenizer.from_pretrained(model_name)
tokenizer.pad_token = tokenizer.eos_token
tokenizer.padding_side = "right"


def preprocess_function(examples):
    inputs = []
    
    for prompt, response in zip(examples["Prompt"], examples["Response"]):
        formatted_text = f"{prompt} [JSON_OUTPUT]\n{response}"
        inputs.append(formatted_text)
    result = tokenizer(inputs, max_length=2048, truncation=True, padding="max_length")
    
    result["labels"] = result["input_ids"].copy()
    
    return result


tokenized_datasets = dataset.map(
    preprocess_function,
    batched=True,
    remove_columns=dataset["train"].column_names,
)

model = AutoModelForCausalLM.from_pretrained(
    model_name,
    device_map="auto",
    torch_dtype=torch.float16,
)


lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["q_proj", "v_proj", "k_proj", "o_proj", "gate_proj", "up_proj", "down_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type=TaskType.CAUSAL_LM,
)


model = prepare_model_for_kbit_training(model)
model = get_peft_model(model, lora_config)
model.print_trainable_parameters()


data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer, 
    mlm=False,  
)


training_args = TrainingArguments(
    output_dir=output_dir,
    evaluation_strategy="epoch",
    learning_rate=2e-4,
    weight_decay=0.01,
    num_train_epochs=3,
    lr_scheduler_type="cosine",
    warmup_ratio=0.1,
    per_device_train_batch_size=1,
    per_device_eval_batch_size=1,
    gradient_accumulation_steps=1,
    save_strategy="epoch",
    save_total_limit=2,
    fp16=True,
    load_best_model_at_end=True,
    metric_for_best_model="loss",
    greater_is_better=False,
    push_to_hub=False,
    report_to="none",
)


trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["test"],
    data_collator=data_collator,
)


trainer.train()


model.save_pretrained(output_dir)
tokenizer.save_pretrained(output_dir)