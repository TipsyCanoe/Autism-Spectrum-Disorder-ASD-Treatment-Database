import "@testing-library/jest-dom";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom"; // Import MemoryRouter
import HomePage from "./HomePage.jsx";

// Mock fetch globally
global.fetch = jest.fn();

describe("HomePage Component", () => {
  beforeEach(() => {
    // Reset fetch mock before each test
    fetch.mockClear();
  });

  test("renders hero section with title and description", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByRole("heading", {
        name: /Welcome to the Autism Resources Database/i,
        level: 1,
      })
    ).toBeInTheDocument();
    expect(
      screen.getByText(
        /A comprehensive database of resources for professionals and families/i
      )
    ).toBeInTheDocument();
    expect(screen.getByText(/Sendan Center/i)).toBeInTheDocument();
  });

  test("renders search section with search button and update database button", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByRole("heading", {
        name: /Find the Resources You Need/i,
        level: 2,
      })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("button", { name: /Search/i })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("button", { name: /Update Database/i })
    ).toBeInTheDocument();
  });

  describe("handleRunJob function", () => {
    test("successfully calls API and displays success message", async () => {
      fetch.mockResolvedValueOnce({
        ok: true,
        json: async () => ({ message: "Job completed successfully", details: "100 entries added" }),
      });

      render(
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      );

      const updateButton = screen.getByRole("button", { name: /Update Database/i });
      fireEvent.click(updateButton);

      expect(updateButton).toBeDisabled();
      expect(screen.getByText(/Updating database.../i)).toBeInTheDocument();

      await waitFor(() => {
        expect(
          screen.getByText(/Job completed successfully Details: 100 entries added/i)
        ).toBeInTheDocument();
      });
      expect(updateButton).not.toBeDisabled();
    });

    test("handles API error and displays error message", async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        json: async () => ({ message: "Job failed", error: "Database connection error" }),
      });

      render(
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      );

      const updateButton = screen.getByRole("button", { name: /Update Database/i });
      fireEvent.click(updateButton);

      expect(updateButton).toBeDisabled();
      expect(screen.getByText(/Updating database.../i)).toBeInTheDocument();

      await waitFor(() => {
        expect(
          screen.getByText(/Error: Job failed Details: Database connection error/i)
        ).toBeInTheDocument();
      });
      expect(updateButton).not.toBeDisabled();
    });

    test("handles network error during API call and displays fallback message", async () => {
      fetch.mockRejectedValueOnce(new Error("Network failure"));

      render(
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      );

      // Mock console.error to verify it's called
      const consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});

      const updateButton = screen.getByRole("button", { name: /Update Database/i });
      fireEvent.click(updateButton);

      expect(updateButton).toBeDisabled();
      expect(screen.getByText(/Updating database.../i)).toBeInTheDocument();

      await waitFor(() => {
        expect(
          screen.getByText(/Failed to trigger job. Check console for details./i)
        ).toBeInTheDocument();
      });
      expect(updateButton).not.toBeDisabled();
      expect(consoleErrorSpy).toHaveBeenCalledWith("Failed to trigger API job:", new Error("Network failure"));
      
      // Restore console.error
      consoleErrorSpy.mockRestore();
    });

    test("displays success message without details if details are not provided", async () => {
      fetch.mockResolvedValueOnce({
        ok: true,
        json: async () => ({ message: "Job completed successfully" }), // No details
      });

      render(
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      );

      const updateButton = screen.getByRole("button", { name: /Update Database/i });
      fireEvent.click(updateButton);

      await waitFor(() => {
        expect(
          screen.getByText(/Job completed successfully/i)
        ).toBeInTheDocument();
         // Ensure "Details:" is not part of the message
        expect(screen.queryByText(/Details:/i)).not.toBeInTheDocument();
      });
    });

    test("displays error message without details if error details are not provided", async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        json: async () => ({ message: "Job failed" }), // No error details
      });

      render(
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      );

      const updateButton = screen.getByRole("button", { name: /Update Database/i });
      fireEvent.click(updateButton);

      await waitFor(() => {
        expect(
          screen.getByText(/Error: Job failed/i)
        ).toBeInTheDocument();
        // Ensure "Details:" is not part of the message
        expect(screen.queryByText(/Details:/i)).not.toBeInTheDocument();
      });
    });
  });

  test("renders features section", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByRole("heading", { name: /Features/i, level: 2 })
    ).toBeInTheDocument();
    expect(screen.getByText(/Extensive Database/i)).toBeInTheDocument();
    expect(screen.getByText(/Professional Guidance/i)).toBeInTheDocument();
    expect(screen.getByText(/Family Support/i)).toBeInTheDocument();
  });

  test("renders testimonials section", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByRole("heading", { name: /What Our Users Say/i, level: 2 })
    ).toBeInTheDocument();
    expect(
      screen.getByText(/"This database has been a lifesaver for our family/i)
    ).toBeInTheDocument();
    expect(
      screen.getByText(/"As a professional, I rely on this database/i)
    ).toBeInTheDocument();
  });

  test("renders other resources section", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByRole("heading", { name: /Other Autism Resources/i, level: 2 })
    ).toBeInTheDocument();
    // Be more specific: look for the heading elements for each resource
    expect(
      screen.getByRole("heading", { name: /Autism Speaks/i, level: 3 })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { name: /National Autism Association/i, level: 3 })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { name: /Autism Society/i, level: 3 })
    ).toBeInTheDocument();
  });

  test("renders footer", () => {
    render(
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    );
    expect(
      screen.getByText(/© 2025 Autism Resources Database. All rights reserved./i)
    ).toBeInTheDocument();
  });
});
