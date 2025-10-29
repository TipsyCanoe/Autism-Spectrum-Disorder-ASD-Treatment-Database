import {
  act,
  fireEvent,
  render,
  screen,
  waitFor,
} from "@testing-library/react";
import { useEffect, useState } from "react";
import useSearch from "./useSearch";

// Mock global fetch
global.fetch = jest.fn();

// Helper Test Component to use the hook and expose its state/functions
const TestComponent = ({
  onStateChange,
  initialQuery = "",
  autoFetchResults = false,
  initialSelectedOptions = null,
  triggerFetchResults = false,
  triggerFetchFilters = false,
}) => {
  const hookResult = useSearch();
  const [internalQuery, setInternalQuery] = useState(initialQuery);

  useEffect(() => {
    if (onStateChange) {
      onStateChange(hookResult);
    }
  }, [hookResult, onStateChange]);

  useEffect(() => {
    if (
      initialSelectedOptions &&
      JSON.stringify(hookResult.selectedOptions) !==
        JSON.stringify(initialSelectedOptions)
    ) {
      act(() => {
        hookResult.setSelectedOptions(initialSelectedOptions);
      });
    }
  }, [
    initialSelectedOptions,
    hookResult.selectedOptions,
    hookResult.setSelectedOptions,
  ]);

  useEffect(() => {
    if (autoFetchResults) {
      act(() => {
        hookResult.fetchResults(internalQuery);
      });
    }
  }, [autoFetchResults, internalQuery, hookResult.fetchResults]);

  useEffect(() => {
    if (triggerFetchResults) {
      act(() => {
        hookResult.fetchResults(internalQuery);
      });
    }
  }, [triggerFetchResults, internalQuery, hookResult.fetchResults]);

  useEffect(() => {
    if (triggerFetchFilters) {
      act(() => {
        hookResult.fetchFilters();
      });
    }
  }, [triggerFetchFilters, hookResult.fetchFilters]);

  return (
    <div>
      <div data-testid="isLoading">{String(hookResult.isLoading)}</div>
      <div data-testid="error">{hookResult.error || "null"}</div>
      <div data-testid="results">{JSON.stringify(hookResult.results)}</div>
      <div data-testid="selectedOptions">
        {JSON.stringify(hookResult.selectedOptions)}
      </div>
      <div data-testid="availableFilters">
        {JSON.stringify(hookResult.availableFilters)}
      </div>
      <input
        type="text"
        data-testid="query-input"
        value={internalQuery}
        onChange={(e) => setInternalQuery(e.target.value)}
      />
      <button onClick={() => hookResult.fetchResults(internalQuery)}>
        Fetch Results
      </button>
      <button
        onClick={() =>
          hookResult.setSelectedOptions(["age:child", "symptom:test"])
        }
      >
        Set Default Test Options
      </button>
      <button onClick={() => hookResult.fetchFilters()}>Refetch Filters</button>
    </div>
  );
};

// Top-level describe block for the hook
describe("useSearch Hook", () => {
  beforeEach(() => {
    fetch.mockClear();
    // Default mock for initial calls on mount
    fetch.mockImplementation((url) => {
      if (url.includes("/api/filters")) {
        return Promise.resolve({
          ok: true,
          json: async () => ({
            age: [],
            symptom: [],
            gender: [],
            medication: [],
          }),
        });
      } else if (url.includes("/api/initial-results")) {
        return Promise.resolve({
          ok: true,
          json: async () => [], // Mock initial results as empty array
        });
      }
      return undefined; // Let specific tests override for other URLs
    });
  });

  test("should initialize with correct default state", async () => {
    render(<TestComponent />);

    // Wait for initial filter loading to complete
    await waitFor(() => {
      expect(screen.getByTestId("availableFilters").textContent).toContain(
        "age"
      );
    });

    // Wait for isLoading to be false after all initial loads are done
    await waitFor(() => {
      expect(screen.getByTestId("isLoading").textContent).toBe("false");
    });

    expect(screen.getByTestId("selectedOptions").textContent).toBe(
      JSON.stringify([])
    );
    expect(screen.getByTestId("results").textContent).toBe(JSON.stringify([]));
    expect(screen.getByTestId("error").textContent).toBe("null");
  });

  describe("fetchFilters", () => {
    test("should fetch and set available filters successfully", async () => {
      const mockFilters = {
        age: ["child", "adult"],
        symptom: ["test"],
        gender: [],
        medication: [],
      };

      fetch.mockImplementation((url) => {
        // Override for this test
        if (url.includes("/api/filters")) {
          return Promise.resolve({
            ok: true,
            json: async () => mockFilters,
          });
        } else if (url.includes("/api/initial-results")) {
          return Promise.resolve({
            ok: true,
            json: async () => [],
          });
        }
        return undefined;
      });

      render(<TestComponent triggerFetchFilters={true} />);

      await waitFor(() => {
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockFilters)
        );
      });
      expect(screen.getByTestId("error").textContent).toBe("null");
    });

    test("should set error state if fetching filters fails", async () => {
      fetch.mockImplementation((url) => {
        // Override for this test
        if (url.includes("/api/filters")) {
          return Promise.resolve({
            ok: false,
            status: 500,
          });
        } else if (url.includes("/api/initial-results")) {
          return Promise.resolve({
            ok: true,
            json: async () => [],
          });
        }
        return undefined;
      });

      const consoleErrorSpy = jest
        .spyOn(console, "error")
        .mockImplementation(() => {});

      render(<TestComponent triggerFetchFilters={true} />);

      await waitFor(() => {
        expect(screen.getByTestId("error").textContent).toBe(
          "Failed to load filters. Please refresh the page."
        );
      });

      // availableFilters might still have the default empty medication array
      expect(screen.getByTestId("availableFilters").textContent).toContain(
        "medication"
      );

      expect(consoleErrorSpy).toHaveBeenCalledWith(
        "Error fetching filters:",
        expect.any(Error)
      );
      consoleErrorSpy.mockRestore();
    });
  });

  describe("fetchResults", () => {
    const mockInitialFiltersForFetchResults = {
      // Renamed to avoid conflict if any
      age: ["adult"],
      symptom: ["anxiety"],
      gender: ["female"],
      medication: ["medication1", "medication2"],
    };

    // This beforeEach is specific to the "fetchResults" describe block
    beforeEach(() => {
      fetch.mockImplementation((url) => {
        if (url.includes("/api/filters")) {
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        } else if (url.includes("/api/initial-results")) {
          return Promise.resolve({
            ok: true,
            json: async () => [],
          });
        }
        // /api/search will be mocked in individual tests within this block
        return undefined;
      });
    });

    test("should fetch and set results successfully with no query or filters", async () => {
      const mockSearchResults = [{ id: 1, title: "Study 1" }];
      fetch.mockImplementation((url) => {
        // Further override for /api/search
        if (url.includes("/api/filters"))
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        if (url.includes("/api/initial-results"))
          return Promise.resolve({ ok: true, json: async () => [] });
        if (url.includes("/api/search"))
          return Promise.resolve({
            ok: true,
            json: async () => mockSearchResults,
          });
        return undefined;
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toContain(
          "adult"
        )
      );
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      await waitFor(() =>
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults)
        )
      );
      expect(screen.getByTestId("isLoading").textContent).toBe("false");
      expect(screen.getByTestId("error").textContent).toBe("null");
    });

    test("should fetch results with a query", async () => {
      const mockSearchResults = [{ id: 2, title: "Query Study" }];
      fetch.mockImplementation((url) => {
        if (url.includes("/api/filters"))
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        if (url.includes("/api/initial-results"))
          return Promise.resolve({ ok: true, json: async () => [] });
        if (url.includes("/api/search"))
          return Promise.resolve({
            ok: true,
            json: async () => mockSearchResults,
          });
        return undefined;
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toContain(
          "adult"
        )
      );
      fireEvent.change(screen.getByTestId("query-input"), {
        target: { value: "test query" },
      });
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      await waitFor(() =>
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults)
        )
      );
    });

    test("should fetch results with selected options", async () => {
      const mockSearchResults = [{ id: 3, title: "Filtered Study" }];
      fetch.mockImplementation((url) => {
        if (url.includes("/api/filters"))
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        if (url.includes("/api/initial-results"))
          return Promise.resolve({ ok: true, json: async () => [] });
        if (url.includes("/api/search"))
          return Promise.resolve({
            ok: true,
            json: async () => mockSearchResults,
          });
        return undefined;
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toContain(
          "adult"
        )
      );
      fireEvent.click(
        screen.getByRole("button", { name: /Set Default Test Options/i })
      );
      await waitFor(() =>
        expect(screen.getByTestId("selectedOptions").textContent).toBe(
          JSON.stringify(["age:child", "symptom:test"])
        )
      );
      fireEvent.change(screen.getByTestId("query-input"), {
        target: { value: "filter query" },
      });
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      await waitFor(() =>
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults)
        )
      );
    });

    test("should set isLoading to true during fetch and false after", async () => {
      fetch.mockImplementation((url) => {
        if (url.includes("/api/filters"))
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        if (url.includes("/api/initial-results"))
          return Promise.resolve({ ok: true, json: async () => [] });
        if (url.includes("/api/search"))
          return new Promise((resolve) =>
            setTimeout(
              () => resolve({ ok: true, json: async () => [] }),
              50
            )
          );
        return undefined;
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toContain(
          "adult"
        )
      );
      act(() => {
        fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      });
      expect(screen.getByTestId("isLoading").textContent).toBe("true");
      await waitFor(
        () => expect(screen.getByTestId("isLoading").textContent).toBe("false"),
        { timeout: 2000 }
      );
    });

    test("should set error state if fetching results fails", async () => {
      fetch.mockImplementation((url) => {
        if (url.includes("/api/filters"))
          return Promise.resolve({
            ok: true,
            json: async () => mockInitialFiltersForFetchResults,
          });
        if (url.includes("/api/initial-results"))
          return Promise.resolve({ ok: true, json: async () => [] });
        if (url.includes("/api/search"))
          return Promise.resolve({ ok: false, status: 500 });
        return undefined;
      });
      const consoleErrorSpy = jest
        .spyOn(console, "error")
        .mockImplementation(() => {});

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toContain(
          "adult"
        )
      );
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      await waitFor(() =>
        expect(screen.getByTestId("error").textContent).toBe(
          "Failed to load results. Please try again."
        )
      );
      expect(screen.getByTestId("results").textContent).toBe(
        JSON.stringify([])
      );
      expect(screen.getByTestId("isLoading").textContent).toBe("false");
      expect(consoleErrorSpy).toHaveBeenCalledWith(
        "Error fetching results:",
        expect.any(Error)
      );
      consoleErrorSpy.mockRestore();
    });
  });

  test("setSelectedOptions should update selectedOptions state", async () => {
    render(<TestComponent />);
    await waitFor(() =>
      expect(screen.getByTestId("availableFilters").textContent).toContain(
        "age"
      )
    );
    const newOptions = ["age:child", "symptom:test"];
    fireEvent.click(
      screen.getByRole("button", { name: /Set Default Test Options/i })
    );
    expect(screen.getByTestId("selectedOptions").textContent).toBe(
      JSON.stringify(newOptions)
    );
  });

  test("fetchResults uses updated selectedOptions correctly", async () => {
    const mockInitialFiltersForThisTest = {
      age: ["adult"],
      symptom: ["anxiety"],
      gender: ["female"],
      medication: ["medication1"],
    };
    fetch.mockImplementation((url) => {
      if (url.includes("/api/filters"))
        return Promise.resolve({
          ok: true,
          json: async () => mockInitialFiltersForThisTest,
        });
      if (url.includes("/api/initial-results"))
        return Promise.resolve({ ok: true, json: async () => [] });
      if (url.includes("/api/search") && !url.includes("filters="))
        return Promise.resolve({
          ok: true,
          json: async () => [{ id: 1, title: "Initial" }],
        });
      if (url.includes("/api/search") && url.includes("filters="))
        return Promise.resolve({
          ok: true,
          json: async () => [{ id: 2, title: "Filtered" }],
        });
      return undefined;
    });

    render(<TestComponent />);
    await waitFor(() =>
      expect(screen.getByTestId("availableFilters").textContent).toContain(
        "adult"
      )
    );
    fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
    await waitFor(() =>
      expect(screen.getByTestId("results").textContent).toContain("Initial")
    );
    fireEvent.click(
      screen.getByRole("button", { name: /Set Default Test Options/i })
    );
    await waitFor(() =>
      expect(screen.getByTestId("selectedOptions").textContent).toBe(
        JSON.stringify(["age:child", "symptom:test"])
      )
    );
    fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
    await waitFor(() =>
      expect(screen.getByTestId("results").textContent).toContain("Filtered")
    );
  });
});
