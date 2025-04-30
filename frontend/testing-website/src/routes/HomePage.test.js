import React from "react";
import { render, screen, within } from "@testing-library/react";
import "@testing-library/jest-dom";
import About from "./About.jsx";
import Contact from "./Contact.jsx"; // Import Contact to allow rendering

// --- Test Suite for About Component (Static Content Only) ---
describe("About Component - Static Content Rendering", () => {
  // Note: We are NOT mocking Contact here, just ensuring About renders without crashing.
  // The Contact component itself won't be visible unless its state is triggered,
  // which we are no longer testing in this file.

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
    expect(
      screen.getByText(
        /We believe that reliable information should be accessible/i
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
    expect(
      screen.getByText(
        /process of articles getting added to our database is thoroughly discussed/i
      )
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
    expect(
      screen.getByText(/maintained by a dedicated team of researchers/i)
    ).toBeInTheDocument();
    // Check for a few team members to confirm rendering
    expect(screen.getByText(/James Harle, MD/i)).toBeInTheDocument();
    expect(screen.getByText(/Shameem Ahmed/i)).toBeInTheDocument();
    expect(screen.getByText(/Richard Jefferson/i)).toBeInTheDocument();
    expect(screen.getByText(/Alex Lo/i)).toBeInTheDocument();
    expect(screen.getByText(/Logan Kalloway/i)).toBeInTheDocument();
    expect(screen.getByText(/Cole Oliva/i)).toBeInTheDocument();
    // Check for roles
    expect(screen.getByText(/Database Architect/i)).toBeInTheDocument();
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
    expect(
      screen.getByText(/committed to maintaining the highest standards/i)
    ).toBeInTheDocument();
    // Check for list items
    expect(
      screen.getByText(/Initial screening by 3rd-party scientific databases/i)
    ).toBeInTheDocument();
    expect(
      screen.getByText(/Review by automated queries/i)
    ).toBeInTheDocument();
    expect(
      screen.getByText(/Double checking when categorizing papers/i)
    ).toBeInTheDocument();
    expect(
      screen.getByText(/Regular reassessment to ensure continued relevance/i)
    ).toBeInTheDocument();
  });

  test('renders "Get Involved" section with button (button presence only)', () => {
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
});
