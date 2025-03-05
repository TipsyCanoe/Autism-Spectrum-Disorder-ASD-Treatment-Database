import React, { useState } from "react";
import "../index.css";
import Contact from "./Contact";

const FAQItem = ({ question, answer, isOpen, onClick }) => {
  return (
    <div className="w-full max-w-2xl border border-gray-300 mb-2 bg-white">
      <button
        className="flex justify-between items-center w-full text-left font-semibold text-black p-4 bg-white focus:outline-none"
        onClick={onClick}
      >
        <span>{question}</span>
        <span className="text-xl">{isOpen ? "−" : "+"}</span>
      </button>
      <div
        className={`overflow-hidden transition-all duration-300 ease-in-out bg-white ${
          isOpen ? "max-h-96" : "max-h-0"
        }`}
      >
        <div className="p-4 text-black border-t border-gray-300">{answer}</div>
      </div>
    </div>
  );
};

const FAQ = () => {
  const [openItem, setOpenItem] = useState(null);
  const [showContact, setShowContact] = useState(false);

  const toggleItem = (index) => {
    setOpenItem(openItem === index ? null : index);
  };

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
        "All resources undergo a thorough review process by autism specialists and experts from the Sendan Center to ensure quality and reliability.",
    },
    {
      question: "Is there a cost to access the database?",
      answer:
        "Basic access to the database is free. Some specialized research papers may require a professional subscription.",
    },
    {
      question: "How can I save resources for later reference?",
      answer:
        "Create a free account to save resources to your personal library and organize them into custom folders.",
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
