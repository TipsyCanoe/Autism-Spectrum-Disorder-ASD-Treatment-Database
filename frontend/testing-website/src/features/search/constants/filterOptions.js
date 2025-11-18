export const ageOptions = [
  { value: "0-5", label: "0-5 years" },
  { value: "6-12", label: "6-12 years" },
  { value: "13-17", label: "13-17 years" },
  { value: "18-25", label: "18-25 years" },
  { value: "26-64", label: "26-64 years" },
  { value: "65+", label: "65+ years" },
];

export const symptomOptions = [
  { value: "irritability", label: "Irritability" },
  { value: "adhd", label: "ADHD symptoms" },
  { value: "hyperactivity", label: "Hyperactivity" },
  { value: "social", label: "Social behaviors" },
  { value: "attention-hyperactivity", label: "Attention and hyperactivity" },
  { value: "asd-severity", label: "ASD severity" },
  {
    value: "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance",
    label:
      "Lethargy/social withdrawal, stereotypy, and hyperactivity/noncompliance",
  },
  { value: "anxiety-reactivity", label: "Anxiety and reactivity" },
];

export const genderOptions = [
  { value: "male", label: "Male" },
  { value: "female", label: "Female" },
  { value: "nonbinary", label: "Non-binary" },
];

export const allFilterOptions = {
  age: ageOptions,
  symptom: symptomOptions,
  gender: genderOptions,
};
