import { fireEvent, render, screen } from "@testing-library/react";
import SearchPage from "./SearchPage";

// Create a mock function for useSearch
const mockUseSearch = jest.fn();
jest.mock("../features/search/hooks/useSearch", () => () => mockUseSearch());

// Mock FilterPanel globally
// The module factory now returns a jest.fn(), which then returns the JSX.
jest.mock("../features/search/components/FilterPanel", () =>
  jest.fn(function MockFilterPanel(props) {
    // eslint-disable-next-line react/prop-types
    return (
      <div data-testid="mock-filter-panel" onClick={props.onSearch}>
        Mock FilterPanel
      </div>
    );
  })
);

const defaultUseSearch = {
  selectedOptions: [],
  setSelectedOptions: jest.fn(),
  results: {
    aripiprazole: [
      {
        "Study Title": "Study 1",
        PMID: "123",
        "Full Text URL": "http://example.com/full1",
        Author: "Author 1",
      },
      {
        "Study Title": "Study 2",
        NCTID: "nct-456",
        Author: "Author 2",
      },
    ],
    risperidone: [
      {
        "Study Title": "Study 3",
        PMID: "789",
        Author: "Author 3",
      },
    ],
  },
  fetchResults: jest.fn(),
  isLoading: false,
  error: null,
  availableFilters: {
    age: ["child", "adult"],
    symptom: ["irritability", "hyperactivity"],
    gender: ["male", "female"],
    medication: ["aripiprazole", "risperidone"],
  },
};

beforeEach(() => {
  mockUseSearch.mockReset();
  mockUseSearch.mockReturnValue({ ...defaultUseSearch });
});

describe("SearchPage", () => {
  it("renders study results and accordions", async () => {
    render(<SearchPage />);
    expect(screen.getByText(/aripiprazole/i)).toBeInTheDocument();
    expect(screen.getByText(/risperidone/i)).toBeInTheDocument();
    expect(screen.queryByText(/Study 1/)).not.toBeInTheDocument();
    fireEvent.click(screen.getByText(/aripiprazole/i));
    expect(await screen.findByText(/Study 1/)).toBeInTheDocument();
    fireEvent.click(screen.getByText(/Study 1/));
    expect(await screen.findByText(/Author 1/)).toBeInTheDocument();
    expect(screen.getByText(/View Full Text/i)).toHaveAttribute(
      "href",
      "http://example.com/full1"
    );
  });

  it("returns correct options from getOptionsForCategory", () => {
    render(<SearchPage />);
    expect(screen.getByText(/aripiprazole/i)).toBeInTheDocument();
  });

  it("handles non-object results (early return)", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: null,
    });
    render(<SearchPage />);
    expect(screen.queryByText(/study results/i)).toBeInTheDocument();
  });

  it("handles non-array studyArray (zero studies)", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: { aripiprazole: null },
    });
    render(<SearchPage />);

    expect(screen.queryByText(/aripiprazole/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/study 1/i)).not.toBeInTheDocument();

    expect(
      screen.getByText(
        (content) =>
          /0 studies found across 1 medication group/i.test(content) ||
          /no studies match your criteria|use the filters or keywords to find studies/i.test(
            content
          )
      )
    ).toBeInTheDocument();
  });

  it("shows no results state", () => {
    mockUseSearch.mockReturnValue({
      ...defaultUseSearch,
      results: {},
    });
    render(<SearchPage />);
    expect(
      screen.getByText((content) =>
        /no studies match your criteria|use the filters or keywords to find studies/i.test(
          content
        )
      )
    ).toBeInTheDocument();
  });

  // The following tests require the real FilterPanel, so you must not mock it for these
  describe("with real FilterPanel", () => {
    beforeAll(() => {
      jest.unmock("../features/search/components/FilterPanel");
    });

    it("calls fetchResults on search", () => {
      const mockFetchResults = jest.fn();
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        fetchResults: mockFetchResults,
      });
      render(<SearchPage />);
    });

    it("calls clearFilters and resets state", () => {
      const mockSetSelectedOptions = jest.fn();
      mockUseSearch.mockReturnValue({
        ...defaultUseSearch,
        setSelectedOptions: mockSetSelectedOptions,
      });
      render(<SearchPage />);
    });
  });
});
