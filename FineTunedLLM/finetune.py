from datasets import load_dataset
from transformers import AutoConfig, AutoTokenizer, AutoModelForCausalLM, Trainer, TrainingArguments, pipeline

# Load dataset
dataset = load_dataset("json", data_files="PromptCompletion.json")

# Load tokenizer
tokenizer = AutoTokenizer.from_pretrained("openai-community/gpt2")
tokenizer.pad_token = tokenizer.eos_token
# Tokenize the dataset
def tokenize_function(examples):
    inputs = examples['prompt']
    outputs = examples['completion']
    
    
    model_inputs = tokenizer(inputs, padding="max_length", truncation=True, max_length=512)
    labels = tokenizer(outputs, padding="max_length", truncation=True, max_length=512)
    

    model_inputs['labels'] = labels['input_ids']
    return model_inputs

tokenized_datasets = dataset.map(tokenize_function, batched=True)

config = AutoConfig.from_pretrained("openai-community/gpt2")

# Load pre-trained model
model = AutoModelForCausalLM.from_pretrained("openai-community/gpt2", config = config)

# Define training arguments
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="no",
    learning_rate=5e-5,
    per_device_train_batch_size=8,
    per_device_eval_batch_size=8,
    num_train_epochs=3,
    weight_decay=0.01,
    save_steps=10_000,
    save_total_limit=2,
    logging_dir="./logs",
)

# Set up Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=None,
)


# Fine-tune the model
trainer.train()

# Save the fine-tuned model
model.save_pretrained("./fine_tuned_model")
tokenizer.save_pretrained("./fine_tuned_model")

# Use the fine-tuned model for text generation
generator = pipeline("text-generation", model=model, tokenizer=tokenizer)
output = generator(r"""{"PMID":18381903,"DOI":"10.1176\/Appi.Ajp.2008.07101560","Title":"Efficacy of adjunctive aripiprazole to either valproate or lithium in bipolar mania patients partially nonresponsive to valproate\/lithium monotherapy: a placebo-controlled study.","Abstract":"The authors evaluated the efficacy and safety of adjunctive aripiprazole in bipolar I patients with mania partially nonresponsive to lithium\/valproate monotherapy. This multicenter, randomized, placebo-controlled study included outpatients experiencing a manic or mixed episode (with or without psychotic features). Patients with partial nonresponse to lithium\/valproate monotherapy (defined as a Young Mania Rating Scale total score >\/=16 at the end of phases 1 and 2, with a decrease of <\/=25% between phases) with target serum concentrations of lithium (0.6-1.0 mmol\/liter) or valproate (50-125 mug\/ml) were randomly assigned in a 2:1 ratio to adjunctive aripiprazole (N=253; 15 or 30 mg\/day) or placebo (N=131) for 6 weeks. Mean improvement from baseline in Young Mania Rating Scale total score at week 6 (primary endpoint) was significantly greater with aripiprazole (-13.3) than with placebo (-10.7). Significant improvements in Young Mania Rating Scale total score with aripiprazole versus placebo occurred from week 1 onward. In addition, the mean improvement in Clinical Global Impression Bipolar Version (CGI-BP) severity of illness (mania) score from baseline to week 6 was significantly greater with aripiprazole (-1.9) than with placebo (-1.6). Discontinuation rates due to adverse events were higher with aripiprazole than with placebo (9% versus 5%, respectively). Akathisia was the most frequently reported extrapyramidal symptom-related adverse event and occurred significantly more frequently among those receiving aripiprazole (18.6%) than among those receiving placebo (5.4%). There were no significant differences between treatments in weight change from baseline to week 6 (+0.55 kg and +0.23 kg for aripiprazole and placebo, respectively; last observation carried forward). Adjunctive aripiprazole therapy showed significant improvements in mania symptoms as early as week 1 and demonstrated a tolerability profile similar to that of monotherapy studies.","Authors":"Vieta Eduard, T'joen Caroline, McQuade Robert D, Carson William H, Marcus Ronald N, Sanchez Raymond, Owen Randall, Nameche Laurence","Journal":"The American journal of psychiatry","Keywords":"Adult, Anticonvulsants, Antimanic Agents, Antipsychotic Agents, Aripiprazole, Dose-Response Relationship, Drug, Double-Blind Method, Drug Therapy, Combination, Female, Humans, Lithium Compounds, Male, Middle Aged, Piperazines, Psychiatric Status Rating Scales, Quinolones, Valproic Acid","URL":"https:\/\/www.ncbi.nlm.nih.gov\/pubmed\/18381903","Affiliations":"Clinical Institute of Neuroscience, Hospital Clinic, University of Barcelona, IDIBAPS, Villarroel 170\/Rossello 140, 08036 Barcelona, Spain. evieta@clinic.ub.es"}""", max_new_tokens=200)
print(output)