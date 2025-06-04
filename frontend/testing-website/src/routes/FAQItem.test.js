import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event"; // Import userEvent

import "@testing-library/jest-dom"; // Ensure jest-dom matchers are available
import FAQItem from "./FAQItem.jsx"; // Adjust import path if needed

describe("FAQItem Component", () => {
  const mockQuestion = "What is testing?";
  const mockAnswer = "Testing is important.";
  const mockOnClick = jest.fn(); // Create a mock function for onClick

  test("calls onClick handler when the button is clicked", async () => {
    const user = userEvent.setup(); // Setup user-event
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={false}
        onClick={mockOnClick} // Pass the mock function
      />
    );

    // Find the button (using the question text is a reliable way)
    const button = screen.getByRole("button", { name: mockQuestion });

    // Simulate a user click
    await user.click(button);

    // Assert that the mock function was called
    expect(mockOnClick).toHaveBeenCalledTimes(1);
  });

  test("renders the answer when open", () => {
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={true} // Set isOpen to true
        onClick={mockOnClick}
      />
    );

    // Check if the answer text is now visible
    expect(screen.getByText(mockAnswer)).toBeVisible();
    // Check the container's max-height style
    const answerDiv = screen
      .getByText(mockAnswer)
      .closest("div.overflow-hidden");
    expect(answerDiv).toHaveClass("max-h-96"); // Or whatever class indicates 'open'
  });

  test('renders the "−" icon when open', () => {
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={true} // Set isOpen to true
        onClick={mockOnClick}
      />
    );

    // Check for the minus icon text
    expect(screen.getByText("−")).toBeInTheDocument();
  });

  test("renders the question", () => {
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={false}
        onClick={mockOnClick}
      />
    );

    // Check if the question text is rendered
    expect(screen.getByText(mockQuestion)).toBeInTheDocument();
  });

  test('renders the "+" icon when closed', () => {
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={false}
        onClick={mockOnClick}
      />
    );

    // Check for the plus icon text
    expect(screen.getByText("+")).toBeInTheDocument();
  });

  test("does not render the answer when closed", () => {
    render(
      <FAQItem
        question={mockQuestion}
        answer={mockAnswer}
        isOpen={false}
        onClick={mockOnClick}
      />
    );

    // Check that the answer text is NOT visible initially
    expect(screen.queryByText(mockAnswer)).not.toBeVisible();
    // Or check the container's max-height style if relying on that for hiding
    const answerDiv = screen
      .getByText(mockAnswer)
      .closest("div.overflow-hidden"); // Find the collapsible div
    expect(answerDiv).toHaveClass("max-h-0");
  });
});
