jest.mock("./FAQItem");
jest.mock("./Contact");

import "@testing-library/jest-dom";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import Contact from "./Contact";
import FAQ from "./FAQ.jsx";
import FAQItem from "./FAQItem";

describe("FAQ Component", () => {
  const user = userEvent.setup(); // Setup user-event for simulating user interactions

  beforeEach(() => {
    // Reset all mocks before each test
    jest.resetAllMocks();

    // Set up FAQItem mock implementation
    FAQItem.mockImplementation(({ question, answer, isOpen, onClick }) => (
      <div data-testid={`faq-item-${question.substring(0, 10)}`}>
        <button onClick={onClick}>{question}</button>
        {isOpen && <span>{answer}</span>}
      </div>
    ));

    // Set up Contact mock implementation
    Contact.mockImplementation(({ isOpen, onClose }) => {
      if (isOpen) {
        return (
          <div data-testid="mock-contact">
            <span>Mock Contact is Open</span>
            <button onClick={onClose}>Close Mock Modal</button>
          </div>
        );
      }
      return null;
    });
  });

  test("renders the main heading", () => {
    render(<FAQ />);
    expect(
      screen.getByRole("heading", {
        name: /Frequently Asked Questions/i,
        level: 1,
      })
    ).toBeInTheDocument();
  });

  test("renders all FAQ items initially closed", () => {
    render(<FAQ />);
    // Check if FAQItem mock was called for each expected question
    const expectedQuestions = [
      /What is the Autism Resources Database/i,
      /Who can use this database/i,
      /How are resources vetted/i,
      /Is there a cost to access/i,
      /How can I save resources/i,
    ];

    expectedQuestions.forEach((questionRegex) => {
      // Find the button rendered by the mock FAQItem
      const questionButton = screen.getByRole("button", {
        name: questionRegex,
      });
      expect(questionButton).toBeInTheDocument();

      // Check that the mock FAQItem component was called with isOpen: false for this question
      expect(FAQItem).toHaveBeenCalledWith(
        expect.objectContaining({
          question: expect.stringMatching(questionRegex),
          isOpen: false, // Initially closed
        }),
        {}
      );
    });

    // Verify no answers are visible initially (since mocks only render answer if isOpen)
    const answerElements = screen.queryAllByText(/.+/); // Get all elements with text
    const faqAnswersText = [
      "The Autism Resources Database is a comprehensive collection",
      "Our database is designed for healthcare professionals",
      "All resources are pulled from established and respected",
      "Basic access to the database is free",
      "We hope to implement a feature to save resources",
    ];
    faqAnswersText.forEach((answerText) => {
      const answerElement = answerElements.find((el) =>
        el.textContent.includes(answerText)
      );
      expect(answerElement).toBeUndefined(); // Answers shouldn't be rendered by the mock
    });
  });

  // For the "closes an open FAQ item when it is clicked again" test
  test("closes an open FAQ item when it is clicked again", async () => {
    render(<FAQ />);
    const firstQuestionRegex = /What is the Autism Resources Database/i;
    const firstQuestionButton = screen.getByRole("button", {
      name: firstQuestionRegex,
    });

    await user.click(firstQuestionButton);

    // Reset the mock to clear previous calls
    FAQItem.mockClear();

    await user.click(firstQuestionButton);

    // Find the specific call for our question instead of using "lastCalledWith"
    const firstQuestionCall = FAQItem.mock.calls.find((call) =>
      firstQuestionRegex.test(call[0].question)
    );

    expect(firstQuestionCall[0].isOpen).toBe(false); // Should be closed again
  });

  // Similarly for the "opens a new item" test
  test("opens a new item and closes the previously open item", async () => {
    render(<FAQ />);
    const firstQuestionRegex = /What is the Autism Resources Database/i;
    const secondQuestionRegex = /Who can use this database/i;
    const firstQuestionButton = screen.getByRole("button", {
      name: firstQuestionRegex,
    });
    const secondQuestionButton = screen.getByRole("button", {
      name: secondQuestionRegex,
    });

    // Action 1: Click to open the first item
    await user.click(firstQuestionButton);

    // Reset the mock to clear previous calls
    FAQItem.mockClear();

    // Action 2: Click to open the second item
    await user.click(secondQuestionButton);

    // Find the specific calls for each question
    const firstQuestionCall = FAQItem.mock.calls.find((call) =>
      firstQuestionRegex.test(call[0].question)
    );
    const secondQuestionCall = FAQItem.mock.calls.find((call) =>
      secondQuestionRegex.test(call[0].question)
    );

    expect(firstQuestionCall[0].isOpen).toBe(false); // First item should be closed
    expect(secondQuestionCall[0].isOpen).toBe(true); // Second item should be open
  });

  // --- Contact Modal Interaction Tests (To be added next) ---

  test('renders the "Contact Us" prompt and button', () => {
    render(<FAQ />);
    expect(
      screen.getByText(/Don't see your question here\?/i)
    ).toBeInTheDocument();
    expect(
      screen.getByRole("button", { name: /Contact Us/i })
    ).toBeInTheDocument();
  });
});
