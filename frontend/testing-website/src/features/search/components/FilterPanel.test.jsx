import "@testing-library/jest-dom";
import { fireEvent, render, screen } from "@testing-library/react";
import FilterPanel from "./FilterPanel";

// Mock FilterSection to isolate FilterPanel logic
jest.mock("./FilterSection", () => {
  const MockFilterSection = ({
    title,
    category,
    options,
    selectedOptions,
    onFilterChange,
    isExpanded,
    onToggle,
  }) => (
    <div data-testid={`filter-section-${category}`}>
      <button onClick={onToggle} data-testid={`toggle-${category}`}>
        {title} {isExpanded ? "Collapse" : "Expand"}
      </button>
      {isExpanded &&
        options.map((opt) => (
          <div key={opt.value}>
            <input
              type="checkbox"
              id={`${category}-${opt.value}`}
              checked={selectedOptions.includes(`${category}:${opt.value}`)}
              onChange={() => onFilterChange(category, opt.value)}
            />
            <label htmlFor={`${category}-${opt.value}`}>{opt.label}</label>
          </div>
        ))}
    </div>
  );
  return MockFilterSection;
});

const mockHandleFilterChange = jest.fn();
const mockSetSearchQuery = jest.fn();
const mockClearFilters = jest.fn();
const mockHandleSearch = jest.fn();

const defaultProps = {
  selectedOptions: [],
  handleFilterChange: mockHandleFilterChange,
  searchQuery: "",
  setSearchQuery: mockSetSearchQuery,
  clearFilters: mockClearFilters,
  handleSearch: mockHandleSearch,
  ageOptions: [{ value: "0-3", label: "0-3 Years" }],
  symptomOptions: [{ value: "social", label: "Social" }],
  genderOptions: [{ value: "male", label: "Male" }],
  medicationOptions: [], // Default to no medication options
};

describe("FilterPanel Component", () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  test("renders correctly with initial props", () => {
    render(<FilterPanel {...defaultProps} />);
    expect(screen.getByText("Filters")).toBeInTheDocument();
    expect(screen.getByLabelText("Keywords")).toBeInTheDocument();
    expect(
      screen.getByRole("button", { name: "Clear All" })
    ).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Search" })).toBeInTheDocument();
    expect(screen.getByTestId("filter-section-age")).toBeInTheDocument();
    expect(screen.getByTestId("filter-section-symptom")).toBeInTheDocument();
    expect(screen.getByTestId("filter-section-gender")).toBeInTheDocument();
  });

  test("calls setSearchQuery when keyword input changes", () => {
    render(<FilterPanel {...defaultProps} />);
    const keywordInput = screen.getByLabelText("Keywords");
    fireEvent.change(keywordInput, { target: { value: "test query" } });
    expect(mockSetSearchQuery).toHaveBeenCalledWith("test query");
  });

  test('calls clearFilters when "Clear All" button is clicked', () => {
    render(<FilterPanel {...defaultProps} />);
    fireEvent.click(screen.getByRole("button", { name: "Clear All" }));
    expect(mockClearFilters).toHaveBeenCalledTimes(1);
  });

  test('calls handleSearch when "Search" button is clicked', () => {
    render(<FilterPanel {...defaultProps} />);
    fireEvent.click(screen.getByRole("button", { name: "Search" }));
    expect(mockHandleSearch).toHaveBeenCalledTimes(1);
  });

  test("toggles filter sections expansion", () => {
    render(<FilterPanel {...defaultProps} />);
    const ageToggle = screen.getByTestId("toggle-age");

    // Initially not expanded (as per mock FilterSection, it would show 'Expand')
    expect(screen.getByText("Age Range Expand")).toBeInTheDocument();
    fireEvent.click(ageToggle);
    // Now expanded (mock FilterSection would show 'Collapse')
    expect(screen.getByText("Age Range Collapse")).toBeInTheDocument();
    fireEvent.click(ageToggle);
    expect(screen.getByText("Age Range Expand")).toBeInTheDocument();
  });

  test("renders medication filter section if medicationOptions are provided", () => {
    const propsWithMedication = {
      ...defaultProps,
      medicationOptions: [{ value: "med1", label: "Medication 1" }],
    };
    render(<FilterPanel {...propsWithMedication} />);
    expect(screen.getByTestId("filter-section-medication")).toBeInTheDocument();
  });

  test("does not render medication filter section if medicationOptions are empty or not provided", () => {
    render(<FilterPanel {...defaultProps} medicationOptions={[]} />);
    expect(
      screen.queryByTestId("filter-section-medication")
    ).not.toBeInTheDocument();

    const propsWithoutMedication = { ...defaultProps };
    delete propsWithoutMedication.medicationOptions;
    render(<FilterPanel {...propsWithoutMedication} />);
    expect(
      screen.queryByTestId("filter-section-medication")
    ).not.toBeInTheDocument();
  });

  test("passes correct props to FilterSection components", () => {
    const selected = ["age:0-3"];
    render(<FilterPanel {...defaultProps} selectedOptions={selected} />);

    // Expand the age filter section first
    const ageToggle = screen.getByTestId("toggle-age");
    fireEvent.click(ageToggle); // Click to expand

    expect(screen.getByTestId("filter-section-age")).toBeInTheDocument();
    const ageCheckbox = screen.getByLabelText("0-3 Years"); // Assuming mock renders this label
    expect(ageCheckbox).toBeChecked();
  });
});
