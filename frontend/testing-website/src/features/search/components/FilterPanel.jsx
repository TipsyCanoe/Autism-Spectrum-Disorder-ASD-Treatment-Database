import { useState } from "react";
import FilterSection from "./FilterSection";

const FilterPanel = ({
  selectedOptions,
  handleFilterChange,
  searchQuery,
  setSearchQuery,
  clearFilters,
  handleSearch,
  ageOptions,
  symptomOptions,
  genderOptions,
  medicationOptions, // Assuming medicationOptions might be passed for a new filter section
}) => {
  const [expandedSections, setExpandedSections] = useState({
    age: false,
    symptoms: false,
    gender: false,
  });

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
          options={ageOptions} // Use prop
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.age}
          onToggle={() => toggleSection("age")}
        />

        <FilterSection
          title="Symptoms & Behaviors"
          category="symptom"
          options={symptomOptions} // Use prop
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.symptoms}
          onToggle={() => toggleSection("symptoms")}
        />

        <FilterSection
          title="Gender"
          category="gender"
          options={genderOptions} // Use prop
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.gender}
          onToggle={() => toggleSection("gender")}
        />

        {/* Example for a new Medication Filter Section, if medicationOptions are provided */}
        {medicationOptions && medicationOptions.length > 0 && (
          <FilterSection
            title="Medication"
            category="medication"
            options={medicationOptions}
            selectedOptions={selectedOptions}
            onFilterChange={handleFilterChange}
            isExpanded={expandedSections.medication} // Add 'medication' to expandedSections state if used
            onToggle={() => toggleSection("medication")} // Add 'medication' to toggleSection logic if used
          />
        )}
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
