import "@testing-library/jest-dom";
import { render, screen } from "@testing-library/react";
import React from "react";
import NotFoundPage from "./NotFoundPage.jsx";

describe("NotFoundPage Component", () => {
  test("renders the 404 message", () => {
    render(<NotFoundPage />);

    // Check if the heading text is rendered
    const headingElement = screen.getByRole("heading", {
      name: /404 page not found/i,
    });
    expect(headingElement).toBeInTheDocument();
  });
});
