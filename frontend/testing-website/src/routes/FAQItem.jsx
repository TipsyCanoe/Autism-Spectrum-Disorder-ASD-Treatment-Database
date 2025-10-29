
const FAQItem = ({ question, answer, isOpen, onClick }) => {
  return (
    <div className="w-full max-w-2xl border border-gray-300 mb-2 bg-white">
      <button
        className="flex justify-between items-center w-full text-left font-semibold body-text text-black p-4 bg-white focus:outline-none"
        onClick={onClick}
        // Add aria-expanded for accessibility
        aria-expanded={isOpen}
        // Optionally add aria-controls if the answer div has an id
      >
        <span>{question}</span>
        <span className="text-xl" aria-hidden="true">
          {isOpen ? "−" : "+"}
        </span>{" "}
        {/* Add aria-hidden */}
      </button>
      <div
        // Optionally add an id here to link with aria-controls
        className={`overflow-hidden transition-all duration-300 ease-in-out bg-white ${
          isOpen ? "max-h-96" : "max-h-0"
        }`}
        // Hide content from accessibility tree when closed
        hidden={!isOpen}
      >
        <div className="p-4 body-text text-black border-t border-gray-300">{answer}</div>
      </div>
    </div>
  );
};

export default FAQItem;
