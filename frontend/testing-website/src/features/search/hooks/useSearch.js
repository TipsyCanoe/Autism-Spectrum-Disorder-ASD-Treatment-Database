import { useCallback, useEffect, useState } from "react";

// Define base API URL to target port 5000 for general API calls
const API_BASE_URL = "http://localhost:5000";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]); // This should be an array of treatment objects
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [availableFilters, setAvailableFilters] = useState({
    age: [],
    symptom: [],
    gender: [],
    medication: [],
  });

  const fetchFilters = useCallback(async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/filters`);
      if (!response.ok) {
        throw new Error("Failed to fetch filters");
      }
      const data = await response.json();
      setAvailableFilters(data);
    } catch (err) {
      console.error("Error fetching filters:", err);
      setError("Failed to load filters. Please refresh the page.");
    }
  }, []);

  useEffect(() => {
    fetchFilters();
  }, [fetchFilters]);

  const fetchResults = useCallback(
    async (query = "") => {
      setIsLoading(true);
      setError(null);
      
      try {
        const params = new URLSearchParams();
        params.append("query", query);
        
        selectedOptions.forEach((filter) => {
          params.append("filters", filter);
        });

        const response = await fetch(
          `${API_BASE_URL}/api/search?${params.toString()}`
        );
        
        if (!response.ok) {
          throw new Error("Failed to fetch search results");
        }
        
        const data = await response.json();

        setResults(data.results || data); 
        
      } catch (err) {
        console.error("Error fetching results:", err);
        setError("Failed to load results. Please try again.");
      } finally {
        setIsLoading(false);
      }
    },
    [selectedOptions]
  );

  const transformApiResponseToTreatmentGroups = (apiResults) => {
    if (Array.isArray(apiResults) && apiResults.length > 0 && apiResults[0].treatment) {
      // Group studies by treatment
      const treatmentMap = new Map();
      
      apiResults.forEach(study => {
        const treatmentName = study.treatment || 'Unknown Treatment';
        
        if (!treatmentMap.has(treatmentName)) {
          treatmentMap.set(treatmentName, {
            treatment: treatmentName,
            studies: []
          });
        }
        
        // Remove treatment field from study object to avoid duplication
        const { treatment, ...studyData } = study;
        treatmentMap.get(treatmentName).studies.push(studyData);
      });
      
      return Array.from(treatmentMap.values());
    }
    
    return apiResults;
  };

  return {
    selectedOptions,
    setSelectedOptions,
    results, // This will be an array of treatment objects: [{ treatment: "...", studies: [...] }]
    fetchResults,
    isLoading,
    error,
    availableFilters,
    fetchFilters,
  };
};

export default useSearch;