from transformers import AutoModelForCausalLM, AutoTokenizer
from peft import PeftModel
import torch
import json
import re

model_path = "./fine_tuned_model"
base_model_path = "meta-llama/Llama-3.2-1B"

tokenizer = AutoTokenizer.from_pretrained(model_path)
base_model = AutoModelForCausalLM.from_pretrained(
    base_model_path, 
    device_map="auto",
    torch_dtype=torch.float16
)
model = PeftModel.from_pretrained(base_model, model_path)
model.eval()

study_text = r"""The authors evaluated the efficacy and safety of adjunctive aripiprazole in bipolar I patients with mania partially nonresponsive to lithium\/valproate monotherapy. This multicenter, randomized, placebo-controlled study included outpatients experiencing a manic or mixed episode (with or without psychotic features). Patients with partial nonresponse to lithium\/valproate monotherapy (defined as a Young Mania Rating Scale total score >\/=16 at the end of phases 1 and 2, with a decrease of <\/=25% between phases) with target serum concentrations of lithium (0.6-1.0 mmol\/liter) or valproate (50-125 mug\/ml) were randomly assigned in a 2:1 ratio to adjunctive aripiprazole (N=253; 15 or 30 mg\/day) or placebo (N=131) for 6 weeks. Mean improvement from baseline in Young Mania Rating Scale total score at week 6 (primary endpoint) was significantly greater with aripiprazole (-13.3) than with placebo (-10.7). Significant improvements in Young Mania Rating Scale total score with aripiprazole versus placebo occurred from week 1 onward. In addition, the mean improvement in Clinical Global Impression Bipolar Version (CGI-BP) severity of illness (mania) score from baseline to week 6 was significantly greater with aripiprazole (-1.9) than with placebo (-1.6). Discontinuation rates due to adverse events were higher with aripiprazole than with placebo (9% versus 5%, respectively). Akathisia was the most frequently reported extrapyramidal symptom-related adverse event and occurred significantly more frequently among those receiving aripiprazole (18.6%) than among those receiving placebo (5.4%). There were no significant differences between treatments in weight change from baseline to week 6 (+0.55 kg and +0.23 kg for aripiprazole and placebo, respectively; last observation carried forward). Adjunctive aripiprazole therapy showed significant improvements in mania symptoms as early as week 1 and demonstrated a tolerability profile similar to that of monotherapy studies."""
prompt = f"give me a json file with the treatment name, study type, study duration, Male:female ratio, age range, treatment dose range, primary outcome area, primary outcome measures, results: primary measure, secondary outcome area, secondary outcome measures, results:secondary measures, tolerability/side effects, safety, drop out rate, and ethnicity percentages. Please add null if you don't know the answer: {study_text} [JSON_OUTPUT]"

inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

outputs = model.generate(
    **inputs,
    max_new_tokens=200,
    temperature=0.7,
    top_p=0.9,
    repetition_penalty=1.2,
    num_beams=4,
    do_sample=True,
    pad_token_id=tokenizer.eos_token_id
)

generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)


raw_json_text = generated_text.split("[JSON_OUTPUT]")[-1].split("[END]")[0].strip()

json_pattern = r'\{.*?\}'
json_matches = re.findall(json_pattern, raw_json_text, re.DOTALL)

if json_matches:
    json_str = json_matches[0]
    try:
        parsed_json = json.loads(json_str)
        print(study_text +"\n\n\n")
        print(json.dumps(parsed_json, indent=2))
    except json.JSONDecodeError:
        print("Found JSON-like text but couldn't parse it")
else:
    print("No JSON found in output")