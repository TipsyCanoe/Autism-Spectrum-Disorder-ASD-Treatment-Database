import { useState, useEffect } from "react";
import FilterPanel from "../features/search/components/FilterPanel";
import useSearch from "../features/search/hooks/useSearch";
import "../index.css";

const SearchPage = () => {
  const {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
    isLoading,
    error,
    availableFilters
  } = useSearch();
  const [searchQuery, setSearchQuery] = useState("");

  useEffect(() => {
    fetchResults("");
  }, []); 

  const handleFilterChange = (category, value) => {
    const filterKey = `${category}:${value}`;
    setSelectedOptions((prev) =>
      prev.includes(filterKey)
        ? prev.filter((item) => item !== filterKey)
        : [...prev, filterKey]
    );
    
    //setTimeout(() => fetchResults(searchQuery), 0);
  };

  const handleSearch = () => {
    fetchResults(searchQuery); 
  };

  const clearFilters = () => {
    setSelectedOptions([]);
    setSearchQuery("");
  };

  const getOptionsForCategory = (category) => {
    if (availableFilters && availableFilters[category] && availableFilters[category].length > 0) {
      return availableFilters[category].map(value => ({
        value,
        label: value.charAt(0).toUpperCase() + value.slice(1).replace(/-/g, ' ')
      }));
    }
    
    switch (category) {
      case 'age': return ageOptions;
      case 'symptom': return symptomOptions;
      case 'gender': return genderOptions;
      case 'medication': return medicationOptions;
      default: return [];
    }
  };

  // Transform the nested JSON structure into a flat array of studies
  const getStudiesFromResults = (results) => {
    if (!results || typeof results !== 'object') return [];
    
    const studies = [];
    Object.entries(results).forEach(([medication, studyArray]) => {
      if (Array.isArray(studyArray)) {
        studyArray.forEach(study => {
          studies.push({
            ...study,
            medication: medication,
            id: `${medication}-${study['Study Title']?.replace(/\s+/g, '-') || Math.random()}`
          });
        });
      }
    });
    return studies;
  };

  const studyList = getStudiesFromResults(results);

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
          handleSearch={handleSearch}
          // Pass the filter options (either from backend or static)
          ageOptions={getOptionsForCategory('age')}
          symptomOptions={getOptionsForCategory('symptom')}
          genderOptions={getOptionsForCategory('gender')}
          medicationOptions={getOptionsForCategory('medication')}
        />

        {/* Results Area */}
        <div className="lg:w-2/3">
          <div className="bg-white p-6 rounded-lg shadow mb-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-bold text-gray-800">Study Results</h2>
              {/* Show count only when not loading, no error, and results exist */}
              {!isLoading && !error && studyList.length > 0 && (
                <p className="text-sm text-gray-500">
                  {studyList.length} studies found
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
                    let label = value;
                    // Find label from options
                    const options = getOptionsForCategory(category);
                    const option = options.find(opt => opt.value === value);
                    if (option) {
                      label = option.label;
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
                {error}
              </p>
            )}

            {/* No Results State (only show if not loading and no error) */}
            {!isLoading && !error && studyList.length === 0 && (
              <div className="text-center py-8">
                <p className="text-gray-500">
                  {selectedOptions.length > 0 || searchQuery
                    ? "No studies match your criteria"
                    : "Use the filters or keywords to find studies"}
                </p>
              </div>
            )}

            {/* Results List (only show if not loading, no error, and results exist) */}
            {!isLoading && !error && studyList.length > 0 && (
              <div className="divide-y divide-gray-200">
                {studyList.map((study) => (
                  <div key={study.id} className="py-6">
                    <div className="flex justify-between items-start mb-2">
                      <h3 className="text-lg font-semibold text-gray-800">
                        {study['Study Title'] || "Untitled Study"}
                      </h3>
                      <span className="inline-block px-3 py-1 text-sm font-medium bg-blue-100 text-blue-800 rounded-full">
                        {study.medication}
                      </span>
                    </div>
                    
                    {/* Study metadata */}
                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4 text-sm text-gray-600">
                      <div>
                        <span className="font-medium">Duration:</span> {study.Duration || 'N/A'}
                      </div>
                      <div>
                        <span className="font-medium">Participants:</span> {study.n || 'N/A'}
                      </div>
                      <div>
                        <span className="font-medium">Age:</span> {study['Age range/mean'] || 'N/A'}
                      </div>
                      <div>
                        <span className="font-medium">M:F Ratio:</span> {study['M:F ratio'] || 'N/A'}
                      </div>
                    </div>

                    {/* Primary outcomes */}
                    <div className="mb-3">
                      <h4 className="font-medium text-gray-800 mb-1">Primary Outcome</h4>
                      <p className="text-sm text-gray-600">
                        <span className="font-medium">{study['Primary Outcome Area']}:</span> {study['Results: Primary measure']}
                      </p>
                      {study['Primary Outcome Measures'] && (
                        <p className="text-xs text-gray-500">Measures: {study['Primary Outcome Measures']}</p>
                      )}
                    </div>

                    {/* Secondary outcomes */}
                    {study['Secondary Outcome Area'] && (
                      <div className="mb-3">
                        <h4 className="font-medium text-gray-800 mb-1">Secondary Outcome</h4>
                        <p className="text-sm text-gray-600">
                          <span className="font-medium">{study['Secondary Outcome Area']}:</span> {study['Results: Secondary Measures']}
                        </p>
                        {study['Secondary Outcome Measures'] && (
                          <p className="text-xs text-gray-500">Measures: {study['Secondary Outcome Measures']}</p>
                        )}
                      </div>
                    )}

                    {/* Treatment details */}
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-3 text-sm">
                      <div>
                        <span className="font-medium text-gray-800">Dose Range:</span>
                        <p className="text-gray-600">{study['Medication/Treatment Dose Range'] || 'N/A'}</p>
                      </div>
                      <div>
                        <span className="font-medium text-gray-800">Drop Out Rate:</span>
                        <p className="text-gray-600">{study['Drop Out Rate'] || 'N/A'}</p>
                      </div>
                    </div>

                    {/* Safety information */}
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
                      <div>
                        <span className="font-medium text-gray-800">Side Effects:</span>
                        <p className="text-gray-600">{study['Tolerability/Side Effects'] || 'N/A'}</p>
                      </div>
                      <div>
                        <span className="font-medium text-gray-800">Safety:</span>
                        <p className="text-gray-600">{study.Safety || 'N/A'}</p>
                      </div>
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
// These will be used as fallbacks if the backend doesn't provide filter options
const ageOptions = [
  { value: "0-5", label: "Infancy/Early Childhood (0-5 years)" },
  { value: "6-12", label: "Childhood (6-12 years)" },
  { value: "13-17", label: "Adolescence (13-17 years)" },
  { value: "18-25", label: "Young Adult (18-25 years)" },
  { value: "26-64", label: "Adult (26-64 years)" },
  { value: "65+", label: "Senior (65+ years)" },
];

const symptomOptions = [
  { value: "irritability", label: "Irritability" },
  { value: "adhd", label: "ADHD symptoms" },
  { value: "hyperactivity", label: "Hyperactivity" },
  { value: "social", label: "Social behaviors" },
  { value: "attention-hyperactivity", label: "Attention and hyperactivity" },
  { value: "asd-severity", label: "ASD severity" },
  { value: "lethargy-withdrawal-stereotypy-hyperactivity-noncompliance", label: "Lethargy/social withdrawal, stereotypy, and hyperactivity/noncompliance" },
  { value: "anxiety-reactivity", label: "Anxiety and reactivity" },
];

const genderOptions = [
  { value: "male", label: "Male" },
  { value: "female", label: "Female" },
  { value: "nonbinary", label: "Non-binary" },
];

const medicationOptions = [
  { value: "aripiprazole", label: "Aripiprazole" },
  { value: "citalopram", label: "Citalopram" },
  { value: "fluoxetine", label: "Fluoxetine" },
  { value: "sertraline", label: "Sertraline" },
  { value: "risperidone", label: "Risperidone" },
];

export default SearchPage;