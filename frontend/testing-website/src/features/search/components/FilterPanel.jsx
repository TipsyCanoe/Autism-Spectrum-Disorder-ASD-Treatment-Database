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
  includeAi,
  setIncludeAi,
}) => {
  const [expandedSections, setExpandedSections] = useState({
    age: false,
    symptoms: false,
    gender: false,
    medication: false, // Add medication to initial state
  });
  
  const [isFilterPanelOpen, setIsFilterPanelOpen] = useState(false);

  const toggleSection = (section) => {
    setExpandedSections((prev) => ({
      ...prev,
      [section]: !prev[section],
    }));
  };

  return (
    <>
      {/* Mobile Filter Toggle Button */}
      <div className="lg:hidden mb-4">
        <button
          onClick={() => setIsFilterPanelOpen(!isFilterPanelOpen)}
          className="w-full flex justify-between items-center px-4 py-3 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
        >
          <span className="font-medium">
            {isFilterPanelOpen ? "Hide Filters" : "Show Filters"}
            {selectedOptions.length > 0 && ` (${selectedOptions.length} active)`}
          </span>
          <svg
            xmlns="http://www.w3.org/2000/svg"
            fill="none"
            viewBox="0 0 24 24"
            strokeWidth={2}
            stroke="currentColor"
            className={`w-5 h-5 transition-transform ${isFilterPanelOpen ? 'rotate-180' : ''}`}
          >
            <path strokeLinecap="round" strokeLinejoin="round" d="M19.5 8.25l-7.5 7.5-7.5-7.5" />
          </svg>
        </button>
      </div>

      {/* Filter Panel */}
      <div className={`
        lg:w-1/3 bg-white p-4 lg:p-6 rounded-lg shadow 
        lg:sticky lg:top-24 lg:self-start
        ${isFilterPanelOpen ? 'block' : 'hidden lg:block'}
      `}>
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-xl font-bold text-gray-800">Filters</h2>
          <button
            onClick={clearFilters}
            className="text-sm text-blue-600 hover:underline"
          >
            Clear All
          </button>
        </div>

        {/* Helpful instructions */}
        <div className="mb-4 p-3 bg-blue-50 border border-blue-200 rounded-lg">
          <p className="text-sm text-gray-700">
            <span className="font-semibold">Tip:</span> Filters and keywords are optional. 
            You can search without selecting any filters. Select multiple options within each category to refine your results.
          </p>
        </div>

        {/* Make the filters area scrollable with fixed height */}
        <div className="overflow-y-auto max-h-[50vh] lg:max-h-[60vh] pr-2">
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
            onKeyDown={(e) => {
              if (e.key === 'Enter') {
                handleSearch();
              }
            }}
          />
        </div>

        {/* AI Results Toggle */}
        <div className="mb-6">
          <label className="flex items-center gap-2 cursor-pointer">
            <input
              type="checkbox"
              checked={includeAi}
              onChange={(e) => setIncludeAi(e.target.checked)}
              className="w-4 h-4 text-navbar-blue rounded border-gray-300 focus:ring-navbar-blue"
            />
            <span className="text-gray-700 font-medium">Include results helped by AI</span>
          </label>
          <p className="text-xs text-gray-500 mt-1 ml-6">
            Uncheck to see only manually verified studies.
          </p>
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

        {/* Gender Filter Section */}
        <FilterSection
          title="Gender"
          category="gender"
          options={genderOptions} // Use prop
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.gender}
          onToggle={() => toggleSection("gender")}
        />

        {/* Symptoms Filter Section */}
        <FilterSection
          title="Symptoms & Behaviors"
          category="symptom"
          options={symptomOptions} // Use prop
          selectedOptions={selectedOptions}
          onFilterChange={handleFilterChange}
          isExpanded={expandedSections.symptoms}
          onToggle={() => toggleSection("symptoms")}
          subtitle="Showing top 15 most common"
        />

        {/* Medication Filter Section */}
        {medicationOptions && medicationOptions.length > 0 && (
          <FilterSection
            title="Medication"
            category="medication"
            options={medicationOptions}
            selectedOptions={selectedOptions}
            onFilterChange={handleFilterChange}
            isExpanded={expandedSections.medication} // Ensure this uses the state
            onToggle={() => toggleSection("medication")} // Ensure this uses the state
            subtitle="Showing top 15 most common"
          />
        )}
      </div>

      {/* Sticky search button always visible */}
      <div className="pt-4 mt-2 border-t border-gray-200 lg:sticky lg:bottom-0 bg-white">
        <button
          onClick={handleSearch}
          className="w-full px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
        >
          Search
        </button>
      </div>
    </div>
    </>
  );
};

export default FilterPanel;
