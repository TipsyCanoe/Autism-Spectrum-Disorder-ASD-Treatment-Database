import { useState } from "react";
import "../../../index.css"; // Adjusted path assuming index.css is at src/index.css
import useSearch from "../hooks/useSearch"; // Adjusted path assuming useSearch is in src/features/search/hooks/
import FilterPanel from "./FilterPanel"; // Corrected path

const SearchPage = () => {
  // Destructure isLoading and error from useSearch
  const {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
    isLoading,
    error,
  } = useSearch();
  const [searchQuery, setSearchQuery] = useState("");

  const handleFilterChange = (category, value) => {
    const filterKey = `${category}:${value}`;
    setSelectedOptions((prev) =>
      prev.includes(filterKey)
        ? prev.filter((item) => item !== filterKey)
        : [...prev, filterKey]
    );
  };

  const handleSearch = () => {
    fetchResults(searchQuery); // Pass the current query
  };

  const clearFilters = () => {
    setSelectedOptions([]);
    setSearchQuery("");
    // Optionally trigger a search with cleared filters/query
    // fetchResults('');
  };

  return (
    <div className="w-11/12 max-w-7xl mx-auto py-8">
      <div className="flex flex-col lg:flex-row gap-8">
        {/* Filters Panel Component */}
        <FilterPanel
          selectedOptions={selectedOptions}
          handleFilterChange={handleFilterChange}
          searchQuery={searchQuery}
          setSearchQuery={setSearchQuery}
          clearFilters={clearFilters}
          handleSearch={handleSearch} // Pass the handler
        />

        {/* Results Area */}
        <div className="lg:w-2/3">
          <div className="bg-white p-6 rounded-lg shadow mb-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-bold text-gray-800">Results</h2>
              {/* Show count only when not loading, no error, and results exist */}
              {!isLoading && !error && results.length > 0 && (
                <p className="text-sm text-gray-500">
                  {results.length} resources found
                </p>
              )}
            </div>

            {/* Active filters display (remains the same) */}
            {selectedOptions.length > 0 && (
              <div className="mb-4">
                <p className="text-sm text-gray-500 mb-2">Active filters:</p>
                <div className="flex flex-wrap gap-2">
                  {selectedOptions.map((filter) => {
                    const [category, value] = filter.split(":");
                    let label = value;
                    // Find label logic (ensure ageOptions etc. are defined below or imported)
                    if (category === "age") {
                      label =
                        ageOptions.find((opt) => opt.value === value)?.label ||
                        value;
                    } else if (category === "symptom") {
                      label =
                        symptomOptions.find((opt) => opt.value === value)
                          ?.label || value;
                    } else if (category === "gender") {
                      label =
                        genderOptions.find((opt) => opt.value === value)
                          ?.label || value;
                    }
                    return (
                      <span
                        key={filter}
                        className="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium bg-gray-200 text-gray-800"
                      >
                        {label}
                        <button
                          onClick={() => handleFilterChange(category, value)}
                          className="ml-1 text-gray-500 hover:text-gray-700"
                          aria-label={`Remove ${label} filter`}
                        >
                          ×
                        </button>
                      </span>
                    );
                  })}
                </div>
              </div>
            )}

            {/* Loading State */}
            {isLoading && (
              <p className="text-center py-8 text-gray-500">Loading...</p>
            )}

            {/* Error State */}
            {error && (
              <p className="text-center py-8 text-red-500">
                Error loading results. Please try again.
              </p>
            )}

            {/* No Results State (only show if not loading and no error) */}
            {!isLoading && !error && results.length === 0 && (
              <div className="text-center py-8">
                <p className="text-gray-500">
                  {selectedOptions.length > 0 || searchQuery
                    ? "No resources match your criteria"
                    : "Use the filters or keywords to find resources"}
                </p>
              </div>
            )}

            {/* Results List (only show if not loading, no error, and results exist) */}
            {!isLoading && !error && results.length > 0 && (
              <div className="space-y-6">
                {results.map((treatmentGroup, groupIndex) => (
                  <div key={groupIndex} className="bg-gray-50 p-4 rounded-lg">
                    <h3 className="text-xl font-bold text-gray-800 mb-4 capitalize">
                      {treatmentGroup.treatment}
                    </h3>
                    <div className="divide-y divide-gray-200 bg-white rounded">
                      {treatmentGroup.studies && treatmentGroup.studies.map((study, studyIndex) => (
                        <div key={studyIndex} className="p-4">
                          <h4 className="text-lg font-semibold text-gray-800 mb-2">
                            {study.title || study["Study Title"] || "Untitled Study"}
                          </h4>
                          <p className="text-sm text-gray-500 mb-2">
                            {study["Publication Date"] && `Published: ${study["Publication Date"]}`}
                            {study.Author && ` • ${study.Author}`}
                            {study.PMID && ` • PMID: ${study.PMID}`}
                          </p>
                          <p className="text-gray-600 mb-3">
                            {study.description || study.Abstract || "No description available"}
                          </p>
                          {study["Primary Outcome Area"] && (
                            <p className="text-sm text-gray-600 mb-2">
                              <span className="font-medium">Outcome:</span> {study["Primary Outcome Area"]}
                            </p>
                          )}
                          {study["Treatment Duration"] && (
                            <p className="text-sm text-gray-600 mb-2">
                              <span className="font-medium">Duration:</span> {study["Treatment Duration"]}
                            </p>
                          )}
                          {study["Similarity Score"] !== undefined && (
                            <p className="text-sm text-gray-500 mb-2">
                              Relevance: {(study["Similarity Score"] * 100).toFixed(0)}%
                            </p>
                          )}
                          {study["Full Text URL"] && (
                            <a
                              href={study["Full Text URL"]}
                              target="_blank"
                              rel="noopener noreferrer"
                              className="inline-block text-blue-600 hover:underline"
                            >
                              View Full Study →
                            </a>
                          )}
                        </div>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

// Define filter options arrays here so they are accessible for the results display logic
// Ensure these are defined if not imported
const ageOptions = [
  { value: "0-3", label: "Early Childhood (0-3 years)" },
  { value: "4-11", label: "School Age (4-11 years)" },
  { value: "12-17", label: "Adolescent (12-17 years)" },
  { value: "18+", label: "Adult (18+ years)" },
];

const symptomOptions = [
  { value: "social", label: "Social Communication" },
  { value: "repetitive", label: "Repetitive Behaviors" },
  { value: "sensory", label: "Sensory Processing" },
  { value: "attention", label: "Attention/Focus" },
  { value: "anxiety", label: "Anxiety" },
];

const genderOptions = [
  { value: "male", label: "Male" },
  { value: "female", label: "Female" },
  { value: "nonbinary", label: "Non-binary" },
];

export default SearchPage;
