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
    <div className="w-11/12 max-w-7xl mx-auto py-8">
      <div className="flex flex-col lg:flex-row gap-8">
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

        <div className="lg:w-2/3">
          <div className="bg-white p-6 rounded-lg shadow mb-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-bold text-gray-800">Study Results</h2>
              {!isLoading && !error && medicationGroups.length > 0 && (
                <p className="text-sm text-gray-500">
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
                <p className="text-sm text-gray-500 mb-2">Active filters:</p>
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
                          className="w-full flex justify-between items-center text-left py-3 px-4 bg-blue-50 hover:bg-blue-100 rounded-md focus:outline-none"
                        >
                          <h3 className="text-lg font-semibold text-blue-700">
                            {group.medicationName} ({group.studyCount}{" "}
                            {group.studyCount === 1 ? "study" : "studies"})
                          </h3>
                          <span className="text-blue-500">
                            {activeMedicationKey === group.medicationName ? (
                              <svg
                                xmlns="http://www.w3.org/2000/svg"
                                fill="none"
                                viewBox="0 0 24 24"
                                strokeWidth={1.5}
                                stroke="currentColor"
                                className="w-6 h-6"
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
                                className="w-6 h-6"
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
                          <div className="pl-4 mt-2 space-y-1">
                            {group.studies.map((study) => (
                              <div key={study.id} className="py-1">
                                <button
                                  onClick={() => toggleStudyAccordion(study.id)}
                                  className="w-full flex justify-between items-center text-left py-2 px-3 bg-gray-50 hover:bg-gray-100 rounded-md focus:outline-none"
                                >
                                  <h4 className="text-md font-medium text-gray-700">
                                    {study.title || "Untitled Study"}
                                  </h4>
                                  <span className="text-gray-400">
                                    {activeStudyId === study.id ? (
                                      <svg
                                        xmlns="http://www.w3.org/2000/svg"
                                        fill="none"
                                        viewBox="0 0 24 24"
                                        strokeWidth={1.5}
                                        stroke="currentColor"
                                        className="w-5 h-5"
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
                                        className="w-5 h-5"
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
                                  <div className="p-4 mt-1 bg-white border border-gray-200 rounded-b-md shadow-sm">
                                    {Object.entries(study).map(
                                      ([key, value]) => {
                                        const EXCLUDED_KEYS = [
                                          "id",
                                          "title",
                                          "Study Title",
                                        ];
                                        if (EXCLUDED_KEYS.includes(key)) {
                                          return null;
                                        }
                                        if (
                                          typeof value === "object" &&
                                          value !== null
                                        ) {
                                          return (
                                            <div key={key} className="mb-3">
                                              <strong className="text-sm text-gray-700 capitalize">
                                                {key.replace(/_/g, " ")}:
                                              </strong>
                                              <pre className="text-sm text-gray-600 bg-gray-50 p-2 rounded overflow-x-auto">
                                                {JSON.stringify(value, null, 2)}
                                              </pre>
                                            </div>
                                          );
                                        }
                                        return (
                                          <div key={key} className="mb-3">
                                            <strong className="text-sm text-gray-700 capitalize">
                                              {key.replace(/_/g, " ")}:
                                            </strong>
                                            <p className="text-sm text-gray-600">
                                              {String(
                                                value !== null &&
                                                  value !== undefined
                                                  ? value
                                                  : "N/A"
                                              )}
                                            </p>
                                          </div>
                                        );
                                      }
                                    )}
                                    {study["Full Text URL"] && (
                                      <a
                                        href={study["Full Text URL"]}
                                        target="_blank"
                                        rel="noopener noreferrer"
                                        className="mt-2 inline-block text-blue-600 hover:underline text-sm"
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
    </div>
  );
};

export default SearchPage;