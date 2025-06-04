import {
  act,
  fireEvent,
  render,
  screen,
  waitFor,
} from "@testing-library/react";
import { useEffect, useState } from "react"; // Import React and useEffect for the TestComponent
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
      // Only set if different to avoid loops if TestComponent re-renders due to parent
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
  }, [autoFetchResults, internalQuery, hookResult.fetchResults]); // hookResult.fetchResults dependency is important

  useEffect(() => {
    if (triggerFetchResults) {
      // Allows tests to trigger fetchResults via prop change
      act(() => {
        hookResult.fetchResults(internalQuery);
      });
    }
  }, [triggerFetchResults, internalQuery, hookResult.fetchResults]);

  useEffect(() => {
    if (triggerFetchFilters) {
      // Allows tests to trigger fetchFilters via prop change
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

describe("useSearch Hook", () => {
  beforeEach(() => {
    fetch.mockClear();
  });

  afterEach(() => {});

  test("should initialize with correct default state", async () => {
    fetch.mockResolvedValueOnce({
      // For initial fetchFilters
      ok: true,
      json: async () => ({ age: [], symptom: [], gender: [] }),
    });

    render(<TestComponent />);

    await waitFor(() => {
      expect(screen.getByTestId("availableFilters").textContent).toBe(
        JSON.stringify({ age: [], symptom: [], gender: [] })
      );
    });

    expect(screen.getByTestId("selectedOptions").textContent).toBe(
      JSON.stringify([])
    );
    expect(screen.getByTestId("results").textContent).toBe(JSON.stringify([]));
    expect(screen.getByTestId("isLoading").textContent).toBe("false");
    expect(screen.getByTestId("error").textContent).toBe("null");
  });

  describe("fetchFilters", () => {
    test("should fetch and set available filters successfully", async () => {
      const mockFilters = {
        age: ["child", "adult"],
        symptom: ["test"],
        gender: [],
      };
      fetch.mockResolvedValueOnce({
        // For initial fetchFilters
        ok: true,
        json: async () => mockFilters,
      });

      render(<TestComponent />);

      await waitFor(() => {
        expect(fetch).toHaveBeenCalledWith("http://localhost:5000/api/filters");
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockFilters)
        );
      });
      expect(screen.getByTestId("error").textContent).toBe("null");
      expect(fetch).toHaveBeenCalledTimes(1);
    });

    test("should set error state if fetching filters fails", async () => {
      fetch.mockResolvedValueOnce({
        // For initial fetchFilters
        ok: false,
        status: 500,
      });
      const consoleErrorSpy = jest
        .spyOn(console, "error")
        .mockImplementation(() => {});

      render(<TestComponent />);

      await waitFor(() => {
        expect(screen.getByTestId("error").textContent).toBe(
          "Failed to load filters. Please refresh the page."
        );
      });
      expect(screen.getByTestId("availableFilters").textContent).toBe(
        JSON.stringify({ age: [], symptom: [], gender: [] })
      );
      expect(consoleErrorSpy).toHaveBeenCalledWith(
        "Error fetching filters:",
        new Error("Failed to fetch filters")
      );
      consoleErrorSpy.mockRestore();
    });
  });

  describe("fetchResults", () => {
    const mockInitialFilters = {
      age: ["adult"],
      symptom: ["anxiety"],
      gender: ["female"],
    };

    beforeEach(() => {
      // Mock the initial successful call to fetchFilters that happens in useEffect
      fetch.mockResolvedValueOnce({
        ok: true,
        json: async () => mockInitialFilters,
      });
    });

    test("should fetch and set results successfully with no query or filters", async () => {
      const mockSearchResults = { results: [{ id: 1, title: "Study 1" }] };
      fetch.mockResolvedValueOnce({
        // For fetchResults
        ok: true,
        json: async () => mockSearchResults,
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockInitialFilters)
        )
      );

      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));

      await waitFor(() => {
        expect(fetch).toHaveBeenLastCalledWith(
          "http://localhost:5000/api/search?query="
        );
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults.results)
        );
      });
      expect(screen.getByTestId("isLoading").textContent).toBe("false");
      expect(screen.getByTestId("error").textContent).toBe("null");
      expect(fetch).toHaveBeenCalledTimes(2); // initial filters + results
    });

    test("should fetch results with a query", async () => {
      const mockSearchResults = { results: [{ id: 2, title: "Query Study" }] };
      fetch.mockResolvedValueOnce({
        // For fetchResults
        ok: true,
        json: async () => mockSearchResults,
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockInitialFilters)
        )
      );

      fireEvent.change(screen.getByTestId("query-input"), {
        target: { value: "test query" },
      });
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));

      await waitFor(() => {
        expect(fetch).toHaveBeenLastCalledWith(
          "http://localhost:5000/api/search?query=test+query"
        );
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults.results)
        );
      });
      expect(fetch).toHaveBeenCalledTimes(2);
    });

    test("should fetch results with selected options", async () => {
      const mockSearchResults = {
        results: [{ id: 3, title: "Filtered Study" }],
      };
      fetch.mockResolvedValueOnce({
        // For fetchResults
        ok: true,
        json: async () => mockSearchResults,
      });

      // Use the TestComponent's button to set options
      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockInitialFilters)
        )
      );

      // Set options using the button
      fireEvent.click(
        screen.getByRole("button", { name: /Set Default Test Options/i })
      );
      await waitFor(() =>
        expect(screen.getByTestId("selectedOptions").textContent).toBe(
          JSON.stringify(["age:child", "symptom:test"])
        )
      );

      // Input query and fetch
      fireEvent.change(screen.getByTestId("query-input"), {
        target: { value: "filter query" },
      });
      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));

      await waitFor(() => {
        expect(fetch).toHaveBeenLastCalledWith(
          "http://localhost:5000/api/search?query=filter+query&filters=age%3Achild&filters=symptom%3Atest"
        );
        expect(screen.getByTestId("results").textContent).toBe(
          JSON.stringify(mockSearchResults.results)
        );
      });
      expect(fetch).toHaveBeenCalledTimes(2);
    });

    test("should set isLoading to true during fetch and false after", async () => {
      fetch.mockImplementationOnce(() => {
        // For initial fetchFilters
        return Promise.resolve({
          ok: true,
          json: async () => mockInitialFilters,
        });
      });
      fetch.mockImplementationOnce(() => {
        // For fetchResults - delayed
        return new Promise((resolve) =>
          setTimeout(
            () =>
              resolve({
                ok: true,
                json: async () => ({ results: [] }),
              }),
            50
          )
        );
      });

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockInitialFilters)
        )
      );

      // Don't await the click, check isLoading immediately after
      act(() => {
        fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
      });

      expect(screen.getByTestId("isLoading").textContent).toBe("true");
      // Wait for the fetch to complete and isLoading to become false
      await waitFor(
        () => expect(screen.getByTestId("isLoading").textContent).toBe("false"),
        { timeout: 2000 }
      );
    });

    test("should set error state if fetching results fails", async () => {
      fetch.mockResolvedValueOnce({
        // For fetchResults
        ok: false,
        status: 500,
      });
      const consoleErrorSpy = jest
        .spyOn(console, "error")
        .mockImplementation(() => {});

      render(<TestComponent />);
      await waitFor(() =>
        expect(screen.getByTestId("availableFilters").textContent).toBe(
          JSON.stringify(mockInitialFilters)
        )
      );

      fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));

      await waitFor(() => {
        expect(screen.getByTestId("error").textContent).toBe(
          "Failed to load results. Please try again."
        );
      });
      expect(screen.getByTestId("results").textContent).toBe(
        JSON.stringify([])
      );
      expect(screen.getByTestId("isLoading").textContent).toBe("false");
      expect(consoleErrorSpy).toHaveBeenCalledWith(
        "Error fetching results:",
        new Error("Failed to fetch search results")
      );
      consoleErrorSpy.mockRestore();
    });
  });

  test("setSelectedOptions should update selectedOptions state", async () => {
    fetch.mockResolvedValueOnce({
      // For initial fetchFilters
      ok: true,
      json: async () => ({ age: [], symptom: [], gender: [] }),
    });
    render(<TestComponent />);
    await waitFor(() =>
      expect(screen.getByTestId("availableFilters").textContent).toBe(
        JSON.stringify({ age: [], symptom: [], gender: [] })
      )
    );

    const newOptions = ["age:child", "symptom:test"]; // These are set by the "Set Default Test Options" button
    fireEvent.click(
      screen.getByRole("button", { name: /Set Default Test Options/i })
    );

    expect(screen.getByTestId("selectedOptions").textContent).toBe(
      JSON.stringify(newOptions)
    );
  });

  test("fetchResults uses updated selectedOptions correctly", async () => {
    const mockInitialFilters = {
      age: ["adult"],
      symptom: ["anxiety"],
      gender: ["female"],
    };
    fetch.mockResolvedValueOnce({
      // Initial fetchFilters
      ok: true,
      json: async () => mockInitialFilters,
    });
    fetch.mockResolvedValueOnce({
      // First fetchResults call (no options set by button yet)
      ok: true,
      json: async () => ({ results: [{ id: 1, title: "Initial" }] }),
    });
    fetch.mockResolvedValueOnce({
      // Second fetchResults call (with options set by button)
      ok: true,
      json: async () => ({ results: [{ id: 2, title: "Filtered" }] }),
    });

    render(<TestComponent />);
    await waitFor(() =>
      expect(screen.getByTestId("availableFilters").textContent).toBe(
        JSON.stringify(mockInitialFilters)
      )
    );

    // First call to fetchResults (no options set via button yet, uses default empty query)
    fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
    await waitFor(() =>
      expect(fetch).toHaveBeenLastCalledWith(
        "http://localhost:5000/api/search?query="
      )
    );
    await waitFor(() =>
      expect(screen.getByTestId("results").textContent).toBe(
        JSON.stringify([{ id: 1, title: "Initial" }])
      )
    );

    // Set options using the button in TestComponent
    fireEvent.click(
      screen.getByRole("button", { name: /Set Default Test Options/i })
    );
    await waitFor(() =>
      expect(screen.getByTestId("selectedOptions").textContent).toBe(
        JSON.stringify(["age:child", "symptom:test"])
      )
    );

    // Second call to fetchResults, should use new options (and default empty query)
    fireEvent.click(screen.getByRole("button", { name: /Fetch Results/i }));
    await waitFor(() =>
      expect(fetch).toHaveBeenLastCalledWith(
        "http://localhost:5000/api/search?query=&filters=age%3Achild&filters=symptom%3Atest"
      )
    );
    await waitFor(() =>
      expect(screen.getByTestId("results").textContent).toBe(
        JSON.stringify([{ id: 2, title: "Filtered" }])
      )
    );

    expect(fetch).toHaveBeenCalledTimes(3); // filters, results1, results2
  });
});
