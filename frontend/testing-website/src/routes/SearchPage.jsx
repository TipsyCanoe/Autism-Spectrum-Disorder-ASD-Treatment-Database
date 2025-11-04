import { useCallback, useMemo, useState } from "react";
import FilterPanel from "../features/search/components/FilterPanel";
import { allFilterOptions } from "../features/search/constants/filterOptions";
import useSearch from "../features/search/hooks/useSearch";
import "../index.css";

const SearchPage = () => {
  const {
    selectedOptions,
    setSelectedOptions,
    results, // This is now an array of objects
    fetchResults,
    isLoading,
    error,
    availableFilters,
  } = useSearch();
  const [searchQuery, setSearchQuery] = useState("");
  const [activeMedicationKey, setActiveMedicationKey] = useState(null);
  const [activeStudyId, setActiveStudyId] = useState(null);

  const handleFilterChange = (category, value) => {
    const filterKey = `${category}:${value}`;
    setSelectedOptions((prev) =>
      prev.includes(filterKey)
        ? prev.filter((item) => item !== filterKey)
        : [...prev, filterKey]
    );
  };

  const handleSearch = () => {
    fetchResults(searchQuery);
  };

  const clearFilters = () => {
    setSelectedOptions([]);
    setSearchQuery("");
  };

  const getOptionsForCategory = (category) => {
    if (
      availableFilters &&
      availableFilters[category] &&
      availableFilters[category].length > 0
    ) {
      return availableFilters[category].map((value) => ({
        value: value,
        label:
          value.charAt(0).toUpperCase() + value.slice(1).replace(/-/g, " "),
      }));
    }
    return allFilterOptions[category] || [];
  };

  const processResultsForMultiLevelAccordion = useCallback((results) => {
    if (!results || !Array.isArray(results)) return [];

    return results.map((treatmentGroup) => {
      const treatmentName = treatmentGroup.treatment || "Unknown Treatment";
      const studyArray = treatmentGroup.studies || [];

      const capitalizedTreatmentName = treatmentName
        .split(' ')
        .map(word => word.charAt(0).toUpperCase() + word.slice(1))
        .join(' ');

      const studies = studyArray.map((rawStudy, index) => {
        const title = rawStudy["Study Title"] || "Untitled Study";
        const pmidOrNctid = rawStudy.PMID || rawStudy.NCTID;
        const id = pmidOrNctid
          ? `${treatmentName}-${title.replace(/\s+/g, "-")}-${pmidOrNctid}`
          : `${treatmentName}-${title.replace(/\s+/g, "-")}-idx${index}`;
        
        // Filter out unwanted fields
        const filteredStudy = Object.fromEntries(
          Object.entries(rawStudy).filter(([key]) => 
            !['Similarity Score', 'Distance', 'Year'].includes(key)
          )
        );

        return {
          ...filteredStudy,
          title: title,
          id: id,
        };
      });

      return {
        medicationName: capitalizedTreatmentName,
        studies,
        studyCount: studies.length,
      };
    });
  }, []);

  const medicationGroups = useMemo(
    () => processResultsForMultiLevelAccordion(results),
    [results, processResultsForMultiLevelAccordion]
  );

  const toggleMedicationAccordion = (medicationKey) => {
    setActiveMedicationKey((prevKey) =>
      prevKey === medicationKey ? null : medicationKey
    );
    setActiveStudyId(null);
  };

  const toggleStudyAccordion = (studyId) => {
    setActiveStudyId((prevId) => (prevId === studyId ? null : studyId));
  };

  return (
    <div className="w-full px-4 lg:w-11/12 max-w-7xl mx-auto py-4 lg:py-8">
      {/* Disclaimer Banner */}
      <div className="bg-yellow-100 border-l-4 border-yellow-500 text-yellow-800 p-3 lg:p-4 mb-4 lg:mb-6 rounded text-sm lg:text-base">
        <strong>Disclaimer:</strong> This website is in development and does not provide medical advice or recommendations. For medical decisions, consult a qualified healthcare professional.
      </div>
      <div className="flex flex-col lg:flex-row gap-4 lg:gap-8">
        <FilterPanel
          selectedOptions={selectedOptions}
          handleFilterChange={handleFilterChange}
          searchQuery={searchQuery}
          setSearchQuery={setSearchQuery}
          clearFilters={clearFilters}
          handleSearch={handleSearch}
          ageOptions={getOptionsForCategory("age")}
          symptomOptions={getOptionsForCategory("symptom")}
          genderOptions={getOptionsForCategory("gender")}
          medicationOptions={getOptionsForCategory("medication")}
        />

        <div className="w-full lg:w-2/3">
          <div className="bg-white p-4 lg:p-6 rounded-lg shadow mb-6">
            <div className="flex flex-col sm:flex-row sm:justify-between sm:items-center mb-4 gap-2">
              <h2 className="text-lg lg:text-xl font-bold text-gray-800">Study Results</h2>
              {!isLoading && !error && medicationGroups.length > 0 && (
                <p className="text-xs lg:text-sm text-gray-500">
                  {medicationGroups.reduce(
                    (acc, group) => acc + group.studyCount,
                    0
                  )}{" "}
                  studies found across {medicationGroups.length} treatment
                  groups
                </p>
              )}
            </div>

            {/* Active filters display */}
            {selectedOptions.length > 0 && (
              <div className="mb-4">
                <p className="text-xs lg:text-sm text-gray-500 mb-2">Active filters:</p>
                <div className="flex flex-wrap gap-2">
                  {selectedOptions.map((filter) => {
                    const [category, value] = filter.split(":");
                    let label = value;
                    if (allFilterOptions[category]) {
                      label =
                        allFilterOptions[category].find(
                          (opt) => opt.value === value
                        )?.label ||
                        value.charAt(0).toUpperCase() +
                          value.slice(1).replace(/-/g, " ");
                    } else {
                      label =
                        value.charAt(0).toUpperCase() +
                        value.slice(1).replace(/-/g, " ");
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

            {isLoading && (
              <p className="text-center py-8 text-gray-500">Loading...</p>
            )}
            {error && <p className="text-center py-8 text-red-500">{error}</p>}

            {!isLoading && !error && medicationGroups.length === 0 && (
              <div className="text-center py-8">
                <p className="text-gray-500">
                  {selectedOptions.length > 0 || searchQuery
                    ? "No studies match your criteria"
                    : "Use the filters or keywords to find studies"}
                </p>
              </div>
            )}

            {!isLoading && !error && medicationGroups.length > 0 && (
              <div className="divide-y divide-gray-200">
                {medicationGroups.map(
                  (group) =>
                    group.studyCount > 0 && (
                      <div key={group.medicationName} className="py-2">
                        <button
                          onClick={() =>
                            toggleMedicationAccordion(group.medicationName)
                          }
                          className="w-full flex justify-between items-center text-left py-2 lg:py-3 px-3 lg:px-4 bg-blue-50 hover:bg-blue-100 rounded-md focus:outline-none"
                        >
                          <h3 className="text-sm lg:text-lg font-semibold text-blue-700 pr-2">
                            {group.medicationName} ({group.studyCount}{" "}
                            {group.studyCount === 1 ? "study" : "studies"})
                          </h3>
                          <span className="text-blue-500 flex-shrink-0">
                            {activeMedicationKey === group.medicationName ? (
                              <svg
                                xmlns="http://www.w3.org/2000/svg"
                                fill="none"
                                viewBox="0 0 24 24"
                                strokeWidth={1.5}
                                stroke="currentColor"
                                className="w-5 h-5 lg:w-6 lg:h-6"
                              >
                                <path
                                  strokeLinecap="round"
                                  strokeLinejoin="round"
                                  d="M19.5 12h-15"
                                />
                              </svg>
                            ) : (
                              <svg
                                xmlns="http://www.w3.org/2000/svg"
                                fill="none"
                                viewBox="0 0 24 24"
                                strokeWidth={1.5}
                                stroke="currentColor"
                                className="w-5 h-5 lg:w-6 lg:h-6"
                              >
                                <path
                                  strokeLinecap="round"
                                  strokeLinejoin="round"
                                  d="M12 4.5v15m7.5-7.5h-15"
                                />
                              </svg>
                            )}
                          </span>
                        </button>

                        {activeMedicationKey === group.medicationName && (
                          <div className="pl-2 lg:pl-4 mt-2 space-y-1">
                            {group.studies.map((study) => (
                              <div key={study.id} className="py-1">
                                <button
                                  onClick={() => toggleStudyAccordion(study.id)}
                                  className="w-full flex justify-between items-center text-left py-2 px-2 lg:px-3 bg-gray-50 hover:bg-gray-100 rounded-md focus:outline-none"
                                >
                                  <h4 className="text-sm lg:text-md font-medium text-gray-700 pr-2">
                                    {study.title || "Untitled Study"}
                                  </h4>
                                  <span className="text-gray-400 flex-shrink-0">
                                    {activeStudyId === study.id ? (
                                      <svg
                                        xmlns="http://www.w3.org/2000/svg"
                                        fill="none"
                                        viewBox="0 0 24 24"
                                        strokeWidth={1.5}
                                        stroke="currentColor"
                                        className="w-4 h-4 lg:w-5 lg:h-5"
                                      >
                                        <path
                                          strokeLinecap="round"
                                          strokeLinejoin="round"
                                          d="M19.5 12h-15"
                                        />
                                      </svg>
                                    ) : (
                                      <svg
                                        xmlns="http://www.w3.org/2000/svg"
                                        fill="none"
                                        viewBox="0 0 24 24"
                                        strokeWidth={1.5}
                                        stroke="currentColor"
                                        className="w-4 h-4 lg:w-5 lg:h-5"
                                      >
                                        <path
                                          strokeLinecap="round"
                                          strokeLinejoin="round"
                                          d="M12 4.5v15m7.5-7.5h-15"
                                        />
                                      </svg>
                                    )}
                                  </span>
                                </button>

                                {activeStudyId === study.id && (
                                  <div className="p-3 lg:p-4 mt-1 bg-white border border-gray-200 rounded-b-md shadow-sm">
                                    {/* Publication Info */}
                                    {(study["Publication Date"] || study.Author || study.PMID) && (
                                      <div className="mb-3 lg:mb-4 text-xs lg:text-sm text-gray-600 break-words">
                                        {study["Publication Date"] && (
                                          <div className="mb-1">Published: {study["Publication Date"]}</div>
                                        )}
                                        {study.Author && (
                                          <div className="mb-1">Author: {study.Author}</div>
                                        )}
                                        {study.PMID && (
                                          <div className="mb-1">PMID: {study.PMID}</div>
                                        )}
                                      </div>
                                    )}
                                    
                                    {/* Abstract/Description - Collapsible */}
                                    {(study.description || study.Abstract) && (
                                      <details className="mb-3 lg:mb-4">
                                        <summary className="cursor-pointer text-xs lg:text-sm text-gray-700 font-bold hover:text-gray-900">
                                          Abstract (Click to expand)
                                        </summary>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words mt-2 pl-4 border-l-2 border-gray-300">
                                          {study.description || study.Abstract}
                                        </p>
                                      </details>
                                    )}
                                    
                                    {/* Study Type */}
                                    {study["Study Type"] && study["Study Type"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Study Type:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Study Type"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Duration */}
                                    {study["Treatment Duration"] && study["Treatment Duration"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Duration:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Treatment Duration"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Sample Size */}
                                    {study["Sample Size"] && study["Sample Size"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Sample Size:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Sample Size"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* M:F Ratio */}
                                    {study["M:F Ratio"] && study["M:F Ratio"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">M:F Ratio:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["M:F Ratio"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Age Range */}
                                    {study["Age Range/Mean"] && study["Age Range/Mean"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Age Range/Mean:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Age Range/Mean"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Medication/Treatment Dose Range */}
                                    {study["Medication/Treatment Dose Range"] && study["Medication/Treatment Dose Range"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Medication/Treatment Dose Range:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Medication/Treatment Dose Range"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Primary Outcome Area */}
                                    {study["Primary Outcome Area"] && study["Primary Outcome Area"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Primary Outcome Area:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Primary Outcome Area"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Primary Outcome Measure */}
                                    {study["Primary Outcome Measure"] && study["Primary Outcome Measure"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Primary Outcome Measures:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Primary Outcome Measure"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Results: Primary measure */}
                                    {study["Results: Primary measure"] && study["Results: Primary measure"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Results: Primary Measure:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Results: Primary measure"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Secondary Outcome Area */}
                                    {study["Secondary Outcome Area"] && study["Secondary Outcome Area"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Secondary Outcome Area:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Secondary Outcome Area"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Secondary Outcome Measures */}
                                    {study["Secondary Outcome Measures"] && study["Secondary Outcome Measures"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Secondary Outcome Measures:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Secondary Outcome Measures"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Tolerability/Side Effects */}
                                    {study["Tolerability/Side Effects"] && study["Tolerability/Side Effects"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Tolerability/Side Effects:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Tolerability/Side Effects"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Safety */}
                                    {study["Safety"] && study["Safety"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Safety:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Safety"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Drop Out Rate */}
                                    {study["Drop Out Rate"] && study["Drop Out Rate"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Drop Out Rate:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Drop Out Rate"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Race/Ethnicity Percentages */}
                                    {study["Race/Ethnicity Percentages"] && study["Race/Ethnicity Percentages"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Race/Ethnicity Percentages:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Race/Ethnicity Percentages"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Notes */}
                                    {study["Notes"] && study["Notes"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Notes:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Notes"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Bias Assessments - Collapsible Section */}
                                    {(study["Sequence Generation (selection bias)"] || 
                                      study["Allocation Concealment (selection bias)"] ||
                                      study["Outcome Assessors Blinding (detection bias)"] ||
                                      study["Clinician and Participant Blinding (performance bias)"] ||
                                      study["Incomplete outcome data (attrition bias)"] ||
                                      study["Selective outcome reporting (reporting bias)"] ||
                                      study["Notes on Biases"]) && (
                                      <details className="mb-2 lg:mb-3">
                                        <summary className="cursor-pointer text-xs lg:text-sm text-gray-700 font-bold hover:text-gray-900">
                                          Bias Assessment (Click to expand)
                                        </summary>
                                        <div className="mt-2 pl-4 border-l-2 border-gray-300">
                                          {study["Sequence Generation (selection bias)"] && study["Sequence Generation (selection bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Sequence Generation (selection bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Sequence Generation (selection bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Allocation Concealment (selection bias)"] && study["Allocation Concealment (selection bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Allocation Concealment (selection bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Allocation Concealment (selection bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Outcome Assessors Blinding (detection bias)"] && study["Outcome Assessors Blinding (detection bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Outcome Assessors Blinding (detection bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Outcome Assessors Blinding (detection bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Clinician and Participant Blinding (performance bias)"] && study["Clinician and Participant Blinding (performance bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Clinician and Participant Blinding (performance bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Clinician and Participant Blinding (performance bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Incomplete outcome data (attrition bias)"] && study["Incomplete outcome data (attrition bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Incomplete outcome data (attrition bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Incomplete outcome data (attrition bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Selective outcome reporting (reporting bias)"] && study["Selective outcome reporting (reporting bias)"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Selective outcome reporting (reporting bias):</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Selective outcome reporting (reporting bias)"]}</p>
                                            </div>
                                          )}
                                          {study["Notes on Biases"] && study["Notes on Biases"] !== "N/A" && (
                                            <div className="mb-2">
                                              <strong className="text-xs lg:text-sm text-gray-700">Notes on Biases:</strong>
                                              <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Notes on Biases"]}</p>
                                            </div>
                                          )}
                                        </div>
                                      </details>
                                    )}
                                    
                                    {/* Treatment Duration */}
                                    {study["Treatment Duration"] && study["Treatment Duration"] !== "N/A" && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Treatment Duration:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600 break-words">{study["Treatment Duration"]}</p>
                                      </div>
                                    )}
                                    
                                    {/* Similarity Score */}
                                    {study["Similarity Score"] !== undefined && (
                                      <div className="mb-2 lg:mb-3">
                                        <strong className="text-xs lg:text-sm text-gray-700">Relevance Score:</strong>
                                        <p className="text-xs lg:text-sm text-gray-600">
                                          {(study["Similarity Score"] * 100).toFixed(0)}%
                                        </p>
                                      </div>
                                    )}
                                    {study["Full Text URL"] && (
                                      <a
                                        href={study["Full Text URL"]}
                                        target="_blank"
                                        rel="noopener noreferrer"
                                        className="mt-2 inline-block text-blue-600 hover:underline text-xs lg:text-sm"
                                      >
                                        View Full Text
                                      </a>
                                    )}
                                  </div>
                                )}
                              </div>
                            ))}
                          </div>
                        )}
                      </div>
                    )
                )}
              </div>
            )}
          </div>
        </div>
      </div>
            {/* Add extra padding and smaller font for mobile */}
      <style>{`
        @media (max-width: 640px) {
          .search-page-banner { font-size: 0.95rem; padding: 0.5rem; }
          .search-page-results { padding: 0.5rem; }
        }
      `}</style>
    </div>
  );
};

export default SearchPage;