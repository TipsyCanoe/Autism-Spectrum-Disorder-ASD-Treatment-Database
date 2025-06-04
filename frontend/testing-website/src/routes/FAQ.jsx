import { useState } from "react";
import "../index.css";
import Contact from "./Contact";
import FAQItem from "./FAQItem";

const FAQ = () => {
  const [openItem, setOpenItem] = useState(null);
  const [showContact, setShowContact] = useState(false);

  const toggleItem = (index) => {
    setOpenItem(openItem === index ? null : index);
  };

  // just list of FAQ items
  const faqItems = [
    {
      question: "What is the Autism Resources Database?",
      answer:
        "The Autism Resources Database is a comprehensive collection of research papers, treatment methodologies, and support resources for autism spectrum disorder.",
    },
    {
      question: "Who can use this database?",
      answer:
        "Our database is designed for healthcare professionals, therapists, educators, researchers, and families affected by autism.",
    },
    {
      question: "How are resources vetted before being added?",
      answer:
        "All resources are pulled from established and respected scientific databases like PubMed to ensure quality and reliability.",
    },
    {
      question: "Is there a cost to access the database?",
      answer:
        "Basic access to the database is free.",
    },
    {
      question: "How can I save resources for later reference?",
      answer:
        "We hope to implement a feature to save resources in the future. Stay tuned!",
    },
  ];

  return (
    <div className="flex flex-col items-center py-8 px-4 bg-white">
      <h1 className="text-3xl font-bold text-black mb-8 text-center">
        Frequently Asked Questions
      </h1>

      <div className="w-full max-w-2xl flex flex-col items-center">
        {faqItems.map((item, index) => (
          <FAQItem
            key={index}
            question={item.question}
            answer={item.answer}
            isOpen={openItem === index}
            onClick={() => toggleItem(index)}
          />
        ))}
      </div>

      <div className="mt-8 text-center">
        <p className="text-black mb-4">Don't see your question here?</p>
        <button
          onClick={() => setShowContact(true)}
          className="px-6 py-2 bg-white text-black border border-black rounded-lg hover:bg-gray-100"
        >
          Contact Us
        </button>
      </div>

      <Contact
        isOpen={showContact}
        onClose={() => setShowContact(false)}
      />
    </div>
  );
};

export default FAQ;