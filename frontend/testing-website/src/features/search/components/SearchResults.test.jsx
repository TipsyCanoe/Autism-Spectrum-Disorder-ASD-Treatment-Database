import "@testing-library/jest-dom";
import { fireEvent, render, screen } from "@testing-library/react";
import SearchPage from "./SearchResults";

// Mock the useSearch hook
const mockFetchResults = jest.fn();
const mockSetSelectedOptions = jest.fn();
const mockUseSearch = jest.fn();

jest.mock("../hooks/useSearch", () => ({
  __esModule: true,
  default: () => mockUseSearch(),
}));

// Mock the FilterPanel component
jest.mock("./FilterPanel", () => {
  const MockFilterPanel = ({
    handleSearch,
    clearFilters,
    searchQuery,
    setSearchQuery,
    handleFilterChange,
    selectedOptions,
  }) => (
    <div data-testid="filter-panel">
      <input
        type="text"
        data-testid="search-query-input"
        value={searchQuery}
        onChange={(e) => setSearchQuery(e.target.value)}
      />
      <button data-testid="search-button" onClick={handleSearch}>
        Search
      </button>
      <button data-testid="clear-filters-button" onClick={clearFilters}>
        Clear Filters
      </button>
      {/* Simulate a filter change */}
      <button
        data-testid="add-age-filter-button"
        onClick={() => handleFilterChange("age", "0-3")}
      >
        Add Age Filter
      </button>
      {/* Display selected options to allow testing removal */}
      {selectedOptions.map((opt) => {
        const [category, value] = opt.split(":");
        return (
          <button
            key={opt}
            data-testid={`remove-${category}-${value}`}
            onClick={() => handleFilterChange(category, value)}
          >
            Remove {category}:{value}
          </button>
        );
      })}
    </div>
  );
  return MockFilterPanel;
});

describe("SearchPage Component (from SearchResults.jsx)", () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockUseSearch.mockReturnValue({
      selectedOptions: [],
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: false,
      error: null,
    });
  });

  test("renders initial state correctly with FilterPanel", () => {
    render(<SearchPage />);
    expect(screen.getByTestId("filter-panel")).toBeInTheDocument();
    expect(screen.getByText(/Results/i)).toBeInTheDocument();
    expect(
      screen.getByText(/Use the filters or keywords to find resources/i)
    ).toBeInTheDocument();
  });

  test("displays loading state", () => {
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: [],
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: true,
      error: null,
    });
    render(<SearchPage />);
    expect(screen.getByText(/Loading.../i)).toBeInTheDocument();
  });

  test("displays error state", () => {
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: [],
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: false,
      error: { message: "Failed to fetch" },
    });
    render(<SearchPage />);
    expect(
      screen.getByText(/Error loading results. Please try again./i)
    ).toBeInTheDocument();
  });

  test('displays "no resources match your criteria" when filters are active and results are empty', () => {
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: ["age:0-3"],
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: false,
      error: null,
    });
    render(<SearchPage />);
    expect(
      screen.getByText(/No resources match your criteria/i)
    ).toBeInTheDocument();
  });

  test("displays results when available", () => {
    const mockResults = [
      {
        treatment: "Treatment A",
        studies: [
          {
            title: "Study 1",
            description: "Description for study 1",
            "Publication Date": "2023",
            Author: "Smith et al.",
            PMID: "12345678"
          }
        ]
      },
      {
        treatment: "Treatment B",
        studies: [
          {
            "Study Title": "Study 2",
            Abstract: "Abstract for study 2",
            "Publication Date": "2022"
          }
        ]
      }
    ];
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: [],
      setSelectedOptions: mockSetSelectedOptions,
      results: mockResults,
      fetchResults: mockFetchResults,
      isLoading: false,
      error: null,
    });
    render(<SearchPage />);
    expect(screen.getByText("Treatment A")).toBeInTheDocument();
    expect(screen.getByText("Study 1")).toBeInTheDocument();
    expect(screen.getByText("Description for study 1")).toBeInTheDocument();
    expect(screen.getByText(/Published: 2023/)).toBeInTheDocument();
    expect(screen.getByText(/Smith et al./)).toBeInTheDocument();
    expect(screen.getByText("Treatment B")).toBeInTheDocument();
    expect(screen.getByText("Study 2")).toBeInTheDocument();
    expect(
      screen.getByText(`${mockResults.length} resources found`)
    ).toBeInTheDocument();
  });

  test("calls fetchResults when search button in FilterPanel is clicked", () => {
    render(<SearchPage />);
    const searchInput = screen.getByTestId("search-query-input");
    fireEvent.change(searchInput, { target: { value: "autism resources" } });
    fireEvent.click(screen.getByTestId("search-button"));
    expect(mockFetchResults).toHaveBeenCalledWith("autism resources");
  });

  test("calls setSelectedOptions and clears query on clearFilters", () => {
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: ["age:0-3"],
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: false,
      error: null,
    });
    render(<SearchPage />);
    const searchInput = screen.getByTestId("search-query-input");
    fireEvent.change(searchInput, { target: { value: "initial query" } });

    fireEvent.click(screen.getByTestId("clear-filters-button"));
    expect(mockSetSelectedOptions).toHaveBeenCalledWith([]);
    expect(searchInput.value).toBe(""); // searchQuery state should be cleared
  });

  test("handles filter changes from FilterPanel (add filter)", () => {
    render(<SearchPage />);
    fireEvent.click(screen.getByTestId("add-age-filter-button")); // Simulates FilterPanel calling handleFilterChange
    expect(mockSetSelectedOptions).toHaveBeenCalled();
  });

  test("displays active filters and allows removal", () => {
    const initialFilters = ["age:0-3", "symptom:social"];
    mockUseSearch.mockReturnValueOnce({
      selectedOptions: initialFilters,
      setSelectedOptions: mockSetSelectedOptions,
      results: [],
      fetchResults: mockFetchResults,
      isLoading: false,
      error: null,
    });
    render(<SearchPage />);

    // Check if active filters are displayed (using the labels defined in SearchResults.jsx)
    expect(screen.getByText("Early Childhood (0-3 years)")).toBeInTheDocument();
    expect(screen.getByText("Social Communication")).toBeInTheDocument();

    // The active filters are rendered by SearchPage directly.
    const removeAgeFilterButton = screen.getByRole("button", {
      name: /Remove Early Childhood \(0-3 years\) filter/i,
    });
    fireEvent.click(removeAgeFilterButton);

    expect(mockSetSelectedOptions).toHaveBeenCalledTimes(1);
  });
});
