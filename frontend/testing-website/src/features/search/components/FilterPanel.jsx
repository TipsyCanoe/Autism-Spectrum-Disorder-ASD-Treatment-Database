import React, { useState } from "react";
import FilterSection from "./FilterSection"; // Assuming FilterSection is in the same directory

// Define filter options centrally
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

const FilterPanel = ({
  selectedOptions,
  handleFilterChange,
  searchQuery,
  setSearchQuery,
  clearFilters,
  handleSearch,
}) => {
  // Manage expansion state locally within the panel
  const [expandedSections, setExpandedSections] = useState({
    age: false,
    symptoms: false,
    gender: false,
  });

  // Toggle section expansion
  const toggleSection = (section) => {
    setExpandedSections((prev) => ({
      ...prev,
      [section]: !prev[section],
    }));
  };

  return (
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

        {/* Age Filter Section */}
        <FilterSection
          title="Age Range"
          category="age"
          options={ageOptions}
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.age}
          onToggle={() => toggleSection("age")}
        />

        {/* Symptoms Filter Section */}
        <FilterSection
          title="Symptoms & Behaviors"
          category="symptom"
          options={symptomOptions}
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.symptoms}
          onToggle={() => toggleSection("symptoms")}
        />

        {/* Gender Filter Section */}
        <FilterSection
          title="Gender"
          category="gender"
          options={genderOptions}
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.gender}
          onToggle={() => toggleSection("gender")}
        />
      </div>

      {/* Sticky search button always visible */}
      <div className="pt-4 mt-2 border-t border-gray-200 sticky bottom-0 bg-white">
        <button
          onClick={handleSearch}
          className="w-full px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
        >
          Search
        </button>
      </div>
    </div>
  );
};

export default FilterPanel;
