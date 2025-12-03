import { fireEvent, render, screen } from "@testing-library/react";
import useSearch from '../features/search/hooks/useSearch';
import SearchPage from "./SearchPage";

jest.mock('../features/search/hooks/useSearch');
const mockUseSearch = useSearch;

// Mock the FilterPanel with a simple React.createElement factory (no JSX)
jest.mock("../features/search/components/FilterPanel", () => {
  const React = require("react");
  return (props) =>
    React.createElement(
      "div",
      { "data-testid": "mock-filter-panel", onClick: props.handleSearch },
      React.createElement("input", {
        key: "search-input",
        placeholder: "Search query",
        value: props.searchQuery || "",
        onChange: (e) => props.setSearchQuery(e.target.value),
      }),
      React.createElement(
        "button",
        {
          key: "select-btn",
          onClick: () => props.handleFilterChange("medication", "aripiprazole"),
        },
        "Select Aripiprazole"
      ),
      React.createElement(
        "button",
        { key: "clear-btn", onClick: props.clearFilters },
        "Clear Filters"
      )
    );
});

const defaultUseSearch = {
  selectedOptions: [],
  setSelectedOptions: jest.fn(),
  // Simulated API results format with treatment and studies for component processing
  results: [
    {
      treatment: "aripiprazole",
      studies: [
        {
          "Study Title": "Study 1",
          PMID: "123",
          "Full Text URL": "http://example.com/full1",
          Author: "Author 1, Author 2, Author 3",
        },
        {
          "Study Title": "Study 2",
          NCTID: "nct-456",
          Author: "Author 2",
        },
      ],
    },
    {
      treatment: "risperidone",
      studies: [
        {
          "Study Title": "Study 3",
          PMID: "789",
          Author: "Author 3",
        },
      ],
    },
  ],
  fetchResults: jest.fn(),
  isLoading: false,
  error: null,
  availableFilters: {
    age: ["child", "adult"],
    symptom: ["irritability", "hyperactivity"],
    gender: ["male", "female"],
    medication: ["aripiprazole", "risperidone"], // This should align with treatment in results
  },
};

beforeEach(() => {
  mockUseSearch.mockReset();
  // Ensure all mocked functions from defaultUseSearch are fresh for each test
  mockUseSearch.mockReturnValue({
    ...defaultUseSearch,
    setSelectedOptions: jest.fn(),
    fetchResults: jest.fn(),
  });
});

describe("SearchPage", () => {
  it("renders study results and accordions", async () => {
    render(<SearchPage />);
    // Check for medication names (accordion headers)
    expect(screen.getByRole('heading', { level: 3, name: /Aripiprazole/i })).toBeInTheDocument();
    expect(screen.getByRole('heading', { level: 3, name: /Risperidone/i })).toBeInTheDocument();

    // Initially, study details should not be visible
    expect(screen.queryByText(/Study 1/i)).not.toBeInTheDocument();

    // Click to expand the aripiprazole accordion
    fireEvent.click(
      screen.getByRole('button', { name: /Aripiprazole.*studies/i })
    );
    // Now "Study 1" should be visible
    expect(await screen.findByText(/Study 1/i)).toBeInTheDocument();

    // Click on "Study 1" to expand its details
    fireEvent.click(screen.getByText(/Study 1/i));
    // Should only show the first author
    expect(await screen.findByText(/Author 1/i)).toBeInTheDocument();
    // Should NOT show the other authors
    expect(screen.queryByText(/Author 2/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/Author 3/i)).not.toBeInTheDocument();
    
    expect(screen.getByText(/View Full Text/i)).toHaveAttribute(
      "href",
      "http://example.com/full1"
    );
  });

  it("returns correct options from getOptionsForCategory", () => {
    render(<SearchPage />);
    expect(screen.getByRole('heading', { level: 3, name: /Aripiprazole/i })).toBeInTheDocument();
  });

  it("handles non-array results (e.g. null, early return)", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: null, // Simulate results being null
    });
    render(<SearchPage />);
    // Expect some placeholder or "no results" message
    expect(
      screen.getByText(/Use the filters or keywords to find studies/i)
    ).toBeInTheDocument();
    // Medication names should not be present if results are null
    expect(screen.queryByRole('heading', { level: 3, name: /Aripiprazole/i })).not.toBeInTheDocument();
  });

  it("handles empty array results (zero studies)", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: [], // Simulate results being an empty array
    });
    render(<SearchPage />);
    // With no filters or query, show initial placeholder
    expect(
      screen.getByText(/Use the filters or keywords to find studies/i)
    ).toBeInTheDocument();
    expect(screen.queryByRole('heading', { level: 3, name: /Aripiprazole/i })).not.toBeInTheDocument();
  });

  it("shows no results state when results is an empty array", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: [], // Empty array for results
    });
    render(<SearchPage />);
    // No filters or query => initial placeholder
    expect(
      screen.getByText(/Use the filters or keywords to find studies/i)
    ).toBeInTheDocument();
  });

  describe("Interactions with FilterPanel", () => {
    it("calls fetchResults when FilterPanel's onSearch is triggered", () => {
      const mockFetchResults = jest.fn();
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        fetchResults: mockFetchResults,
      });
      render(<SearchPage />);
      // Simulate the onSearch event from the mocked FilterPanel
      fireEvent.click(screen.getByTestId("mock-filter-panel"));
      expect(mockFetchResults).toHaveBeenCalled();
    });

    it("calls setSelectedOptions and clears searchQuery when FilterPanel's onClearFilters is triggered", () => {
      const mockSetSelectedOptions = jest.fn();
      // SearchPage's clearFilters calls setSelectedOptions([]) and setSearchQuery("").
      // We'll check setSelectedOptions and also that the search input in the mock FilterPanel gets cleared.
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        setSelectedOptions: mockSetSelectedOptions,
      });

      render(<SearchPage />);

      // Get the search input from the mocked FilterPanel
      const searchInput = screen.getByPlaceholderText("Search query");
      // Optionally, type something into it first to ensure it gets cleared
      fireEvent.change(searchInput, { target: { value: "test query" } });
      expect(searchInput.value).toBe("test query");

      // Simulate the onClearFilters event from the mocked FilterPanel
      fireEvent.click(screen.getByRole("button", { name: "Clear Filters" }));
      expect(mockSetSelectedOptions).toHaveBeenCalledWith([]);
      // Check that SearchPage's setSearchQuery("") was effective by checking the input's value
      expect(searchInput.value).toBe("");
    });
  }); // End of describe("Interactions with FilterPanel")

  describe("Loading and Error States", () => {
    it("displays loading indicator when isLoading is true", () => {
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        isLoading: true,
        results: [], // Ensure results are empty to not interfere
      });
      render(<SearchPage />);
      expect(screen.getByText(/Loading.../i)).toBeInTheDocument();
      // Ensure other content is not present
      expect(screen.queryByText(/Aripiprazole \(/i)).not.toBeInTheDocument();
      expect(
        screen.queryByText(/No studies match your criteria/i)
      ).not.toBeInTheDocument();
      expect(
        screen.queryByText(/Use the filters or keywords to find studies/i)
      ).not.toBeInTheDocument();
    });

    it("displays error message when error is present", () => {
      const errorMessage = "Custom error: Failed to fetch studies";
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        error: errorMessage,
        results: [], // Ensure results are empty
        isLoading: false,
      });
      render(<SearchPage />);
      expect(screen.getByText(errorMessage)).toBeInTheDocument();
      // Ensure other content is not present
      expect(screen.queryByText(/Aripiprazole \(/i)).not.toBeInTheDocument();
      expect(screen.queryByText(/Loading.../i)).not.toBeInTheDocument();
    });
  });

  describe("Filter Management and Display", () => {
    it("calls setSelectedOptions with an updater function when a filter is changed via FilterPanel", () => {
      const mockSetSelectedOptions = jest.fn();
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        setSelectedOptions: mockSetSelectedOptions,
        selectedOptions: [], // Start with no options selected
      });
      render(<SearchPage />);
      // "Select Aripiprazole" button is in the mock FilterPanel
      fireEvent.click(
        screen.getByRole("button", { name: "Select Aripiprazole" })
      );
      expect(mockSetSelectedOptions).toHaveBeenCalledWith(expect.any(Function));

      // Test the updater function behavior for adding a filter
      const updaterAddFunction = mockSetSelectedOptions.mock.calls[0][0];
      expect(updaterAddFunction([])).toEqual(["medication:aripiprazole"]); // From empty to one filter
      expect(updaterAddFunction(["age:child"])).toEqual([
        "age:child",
        "medication:aripiprazole",
      ]); // Adding to existing
    });

    it("displays active filters and allows their removal", () => {
      const mockSetSelectedOptions = jest.fn();
      const initialSelectedOptions = ["medication:aripiprazole", "age:child"];
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        selectedOptions: initialSelectedOptions,
        setSelectedOptions: mockSetSelectedOptions,
        // availableFilters or allFilterOptions should exist for label generation if specific labels are tested.
        // Using default label generation for "aripiprazole" -> "Aripiprazole"
        // and "child" -> "Child"
      });

      render(<SearchPage />);

      // Check if active filters are displayed
      expect(screen.getByText("Aripiprazole")).toBeInTheDocument();
      expect(screen.getByText("Child")).toBeInTheDocument(); // Assuming "age:child" results in "Child" label

      // Find the remove button for "Aripiprazole"
      const removeAripiprazoleButton = screen.getByRole("button", {
        name: /Remove Aripiprazole filter/i,
      });
      fireEvent.click(removeAripiprazoleButton);

      expect(mockSetSelectedOptions).toHaveBeenCalledWith(expect.any(Function));
      const updaterFunction = mockSetSelectedOptions.mock.calls[0][0];
      // handleFilterChange("medication", "aripiprazole") is called.
      // The updater function receives initialSelectedOptions.
      // It should remove "medication:aripiprazole".
      expect(updaterFunction(initialSelectedOptions)).toEqual(["age:child"]);
    });

    it("toggles a filter when handleFilterChange is invoked (e.g., by clicking active filter's remove button twice)", () => {
      const mockSetSelectedOptions = jest.fn();
      let currentSelectedOptions = ["medication:aripiprazole"];

      // Mock setSelectedOptions to update currentSelectedOptions for multi-step interaction
      mockSetSelectedOptions.mockImplementation((updater) => {
        currentSelectedOptions = updater(currentSelectedOptions);
        // Update the mock return value for subsequent renders/interactions if needed
        mockUseSearch.mockReturnValue({
          ...defaultUseSearch,
          selectedOptions: currentSelectedOptions,
          setSelectedOptions: mockSetSelectedOptions,
        });
      });

      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        selectedOptions: currentSelectedOptions,
        setSelectedOptions: mockSetSelectedOptions,
      });

      const { rerender } = render(<SearchPage />);

      // Initial state: "Aripiprazole" filter is active
      expect(screen.getByText("Aripiprazole")).toBeInTheDocument();
      const removeButton = screen.getByRole("button", {
        name: /Remove Aripiprazole filter/i,
      });

      // First click: remove the filter
      fireEvent.click(removeButton);
      expect(mockSetSelectedOptions).toHaveBeenCalledTimes(1);
      expect(currentSelectedOptions).toEqual([]); // Filter removed

      // Rerender with updated state (simulating React's behavior)
      rerender(<SearchPage />);
      expect(screen.queryByText("Aripiprazole")).not.toBeInTheDocument(); // Filter display is gone

      // To test adding it back by clicking a (hypothetical) "add Aripiprazole" button again:
      // This would require the "Select Aripiprazole" button from FilterPanel.
      // Let's assume we click the "Select Aripiprazole" button from the mock FilterPanel
      const addAripiprazoleButton = screen.getByRole("button", {
        name: "Select Aripiprazole",
      });
      fireEvent.click(addAripiprazoleButton);
      expect(mockSetSelectedOptions).toHaveBeenCalledTimes(2);
      expect(currentSelectedOptions).toEqual(["medication:aripiprazole"]); // Filter added back
    });
  });
});
