import React, { useState } from "react";
import { Outlet } from "react-router-dom";
import useSearch from "../features/search/hooks/useSearch";
import "../index.css";

const SearchPage = () => {
  const { selectedOptions, setSelectedOptions, results, fetchResults } =
    useSearch();
  const [searchQuery, setSearchQuery] = useState("");

  // Filter categories - set all to collapsed by default
  const [expandedSections, setExpandedSections] = useState({
    age: false,
    symptoms: false,
    gender: false,
  });

  // Toggle section expansion
  const toggleSection = (section) => {
    setExpandedSections({
      ...expandedSections,
      [section]: !expandedSections[section],
    });
  };

  // Handle checkbox changes
  const handleFilterChange = (category, value) => {
    const filterKey = `${category}:${value}`;

    if (selectedOptions.includes(filterKey)) {
      setSelectedOptions(selectedOptions.filter((item) => item !== filterKey));
    } else {
      setSelectedOptions([...selectedOptions, filterKey]);
    }
  };

  // Handle search submission
  const handleSearch = () => {
    fetchResults();
  };

  // Clear all filters
  const clearFilters = () => {
    setSelectedOptions([]);
    setSearchQuery("");
  };

  return (
    <div className="w-11/12 max-w-7xl mx-auto py-8">
      <div className="flex flex-col lg:flex-row gap-8">
        {/* Filters Panel */}
        <div className="lg:w-1/3 bg-white p-6 rounded-lg shadow relative">
          <div className="flex justify-between items-center mb-4">
            <h2 className="text-xl font-bold text-gray-800">Filters</h2>
            <button
              onClick={clearFilters}
              className="text-sm text-blue-600 hover:underline"
            >
              Clear All
            </button>
          </div>

          {/* Make the filters area scrollable with fixed height */}
          <div className="overflow-y-auto max-h-[60vh] pr-2">
            {/* Search Query */}
            <div className="mb-6">
              <label
                htmlFor="search-query"
                className="block font-medium text-gray-700 mb-2"
              >
                Keywords
              </label>
              <input
                type="text"
                id="search-query"
                className="w-full p-2 border border-gray-300 rounded"
                placeholder="Search terms..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
              />
            </div>

            {/* Age Filter */}
            <div className="mb-6 border-t border-gray-200 pt-4">
              <button
                className="flex justify-between items-center w-full text-left font-medium text-gray-700 mb-2"
                onClick={() => toggleSection("age")}
              >
                <span>Age Range</span>
                <span>{expandedSections.age ? "−" : "+"}</span>
              </button>

              {expandedSections.age && (
                <div className="ml-2 space-y-2">
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("age:0-3")}
                      onChange={() => handleFilterChange("age", "0-3")}
                    />
                    <span>Early Childhood (0-3 years)</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("age:4-11")}
                      onChange={() => handleFilterChange("age", "4-11")}
                    />
                    <span>School Age (4-11 years)</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("age:12-17")}
                      onChange={() => handleFilterChange("age", "12-17")}
                    />
                    <span>Adolescent (12-17 years)</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("age:18+")}
                      onChange={() => handleFilterChange("age", "18+")}
                    />
                    <span>Adult (18+ years)</span>
                  </label>
                </div>
              )}
            </div>

            {/* Symptoms Filter */}
            <div className="mb-6 border-t border-gray-200 pt-4">
              <button
                className="flex justify-between items-center w-full text-left font-medium text-gray-700 mb-2"
                onClick={() => toggleSection("symptoms")}
              >
                <span>Symptoms & Behaviors</span>
                <span>{expandedSections.symptoms ? "−" : "+"}</span>
              </button>

              {expandedSections.symptoms && (
                <div className="ml-2 space-y-2">
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("symptom:social")}
                      onChange={() => handleFilterChange("symptom", "social")}
                    />
                    <span>Social Communication</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("symptom:repetitive")}
                      onChange={() =>
                        handleFilterChange("symptom", "repetitive")
                      }
                    />
                    <span>Repetitive Behaviors</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("symptom:sensory")}
                      onChange={() => handleFilterChange("symptom", "sensory")}
                    />
                    <span>Sensory Processing</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("symptom:attention")}
                      onChange={() =>
                        handleFilterChange("symptom", "attention")
                      }
                    />
                    <span>Attention/Focus</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("symptom:anxiety")}
                      onChange={() => handleFilterChange("symptom", "anxiety")}
                    />
                    <span>Anxiety</span>
                  </label>
                </div>
              )}
            </div>

            {/* Gender Filter */}
            <div className="mb-6 border-t border-gray-200 pt-4">
              <button
                className="flex justify-between items-center w-full text-left font-medium text-gray-700 mb-2"
                onClick={() => toggleSection("gender")}
              >
                <span>Gender</span>
                <span>{expandedSections.gender ? "−" : "+"}</span>
              </button>

              {expandedSections.gender && (
                <div className="ml-2 space-y-2">
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("gender:male")}
                      onChange={() => handleFilterChange("gender", "male")}
                    />
                    <span>Male</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("gender:female")}
                      onChange={() => handleFilterChange("gender", "female")}
                    />
                    <span>Female</span>
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      className="mr-2"
                      checked={selectedOptions.includes("gender:nonbinary")}
                      onChange={() => handleFilterChange("gender", "nonbinary")}
                    />
                    <span>Non-binary</span>
                  </label>
                </div>
              )}
            </div>
          </div>

          {/* Sticky search button always visible */}
          <div className="pt-4 mt-2 border-t border-gray-200 sticky bottom-0 bg-white">
            <button
              onClick={handleSearch}
              className="w-full px-6 py-2 bg-custom text-white rounded-lg hover:bg-custom-hover"
            >
              Search
            </button>
          </div>
        </div>

        {/* Results Area */}
        <div className="lg:w-2/3">
          <div className="bg-white p-6 rounded-lg shadow mb-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-bold text-gray-800">Results</h2>
              {results.length > 0 && (
                <p className="text-sm text-gray-500">
                  {results.length} resources found
                </p>
              )}
            </div>

            {/* Active filters display */}
            {selectedOptions.length > 0 && (
              <div className="mb-4">
                <p className="text-sm text-gray-500 mb-2">Active filters:</p>
                <div className="flex flex-wrap gap-2">
                  {selectedOptions.map((filter) => {
                    const [category, value] = filter.split(":");
                    return (
                      <span
                        key={filter}
                        className="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium bg-gray-200 text-gray-800"
                      >
                        {category === "age" && `Age: ${value}`}
                        {category === "symptom" && `Symptom: ${value}`}
                        {category === "gender" && `Gender: ${value}`}
                        <button
                          onClick={() => handleFilterChange(category, value)}
                          className="ml-1 text-gray-500 hover:text-gray-700"
                        >
                          ×
                        </button>
                      </span>
                    );
                  })}
                </div>
              </div>
            )}

            {/* Results list remains unchanged */}
            {results.length === 0 ? (
              <div className="text-center py-8">
                <p className="text-gray-500">
                  {selectedOptions.length > 0
                    ? "No resources match your filter criteria"
                    : "Use the filters to find resources"}
                </p>
              </div>
            ) : (
              <div className="divide-y divide-gray-200">
                {results.map((result, index) => (
                  <div key={index} className="py-4">
                    <h3 className="text-lg font-semibold text-gray-800">
                      {result.title || "Untitled Resource"}
                    </h3>
                    <p className="text-sm text-gray-500 mb-2">
                      {result.type || "Resource"}
                      {result.publishDate && ` • ${result.publishDate}`}
                      {result.author && ` • ${result.author}`}
                    </p>
                    <p className="text-gray-600 mb-2">{result.description}</p>
                    <div className="flex gap-2">
                      {result.tags &&
                        result.tags.map((tag) => (
                          <span
                            key={tag}
                            className="inline-block px-2 py-1 text-xs font-medium bg-gray-100 text-gray-700 rounded"
                          >
                            {tag}
                          </span>
                        ))}
                    </div>
                    <a
                      href={result.url}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="mt-2 inline-block text-blue-600 hover:underline"
                    >
                      View Resource
                    </a>
                  </div>
                ))}
              </div>
            )}
          </div>

          {/* Outlet for child routes */}
          <Outlet />
        </div>
      </div>
    </div>
  );
};

export default SearchPage;
