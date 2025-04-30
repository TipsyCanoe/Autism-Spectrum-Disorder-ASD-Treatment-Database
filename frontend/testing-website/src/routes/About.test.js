import React from "react";
import { render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import "@testing-library/jest-dom";
import About from "./About.jsx";
// We still need Contact imported so About can render without error, but we don't mock it.
import Contact from "./Contact.jsx";

// --- Test Suite for About Component (Static Content Only) ---
describe("About Component - Static Content Rendering", () => {
  // --- Static Rendering Tests ---
  const user = userEvent.setup();

  test("renders main heading and mission statement", () => {
    render(<About />);
    expect(
      screen.getByRole("heading", {
        name: /About the Autism Resources Database/i,
        level: 1,
      })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { name: /Our Mission/i, level: 2 })
    ).toBeInTheDocument();
    expect(
      screen.getByText(
        /Our mission is to bridge the gap between research and practice/i
      )
    ).toBeInTheDocument();
  });

  test("renders partnership section", () => {
    render(<About />);
    expect(
      screen.getByRole("heading", {
        name: /Our Partnership with the Sendan Center/i,
        level: 2,
      })
    ).toBeInTheDocument();
    expect(
      screen.getByText(/developed in collaboration with Jim Harle/i)
    ).toBeInTheDocument();
  });

  test('renders "What You\'ll Find" section', () => {
    render(<About />);
    expect(
      screen.getByRole("heading", {
        name: /What You'll Find in Our Database/i,
        level: 2,
      })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { name: /Research Papers/i, level: 3 })
    ).toBeInTheDocument();
    expect(
      screen.getByText(/Established studies on autism interventions/i)
    ).toBeInTheDocument();
  });

  test('renders "Our Team" section', () => {
    render(<About />);
    expect(
      screen.getByRole("heading", { name: /Our Team/i, level: 2 })
    ).toBeInTheDocument();
    // Check for a few team members to confirm rendering
    expect(screen.getByText(/James Harle, MD/i)).toBeInTheDocument();
    expect(screen.getByText(/Cole Oliva/i)).toBeInTheDocument();
    expect(
      screen.getByText(/Frontend and Endpoints Architect/i)
    ).toBeInTheDocument();
  });

  test('renders "Resource Vetting Process" section', () => {
    render(<About />);
    expect(
      screen.getByRole("heading", {
        name: /Our Resource Vetting Process/i,
        level: 2,
      })
    ).toBeInTheDocument();
    // Check for one of the list items
    expect(
      screen.getByText(/Initial screening by 3rd-party scientific databases/i)
    ).toBeInTheDocument();
  });

  test('renders "Get Involved" section with button (presence only)', () => {
    render(<About />);
    const getInvolvedSection = screen
      .getByRole("heading", { name: /Get Involved/i, level: 2 })
      .closest("section");
    expect(getInvolvedSection).toBeInTheDocument();
    expect(
      within(getInvolvedSection).getByText(
        /We welcome contributions from researchers/i
      )
    ).toBeInTheDocument();
    // Only check that the button exists, not its functionality
    expect(
      within(getInvolvedSection).getByRole("button", { name: /Contact Us/i })
    ).toBeInTheDocument();
  });

  test("opens the Contact modal when Contact Us button is clicked", async () => {
    render(<About />);
    const contactButton = screen.getByRole("button", { name: /Contact Us/i });
    await user.click(contactButton);

    // Contact modal should be open now - use a more specific selector
    expect(
      screen.getByRole("heading", { name: /Contact Us/i, level: 2 })
    ).toBeInTheDocument();
    expect(screen.getByLabelText(/Name/i)).toBeInTheDocument();

    // Test closing
    const closeButton = screen.getByRole("button", { name: /Cancel/i });
    await user.click(closeButton);

    // Modal should be gone
    expect(screen.queryByLabelText(/Name/i)).not.toBeInTheDocument();
  });
});
