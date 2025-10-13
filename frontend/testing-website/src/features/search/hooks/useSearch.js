import { useCallback, useEffect, useState } from "react";

// Get API URL from environment variable with fallback
const API_BASE_URL = process.env.REACT_APP_PYTHON_API_URL || "";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]); // This should be an array of treatment objects
  const [isLoading, setIsLoading] = useState(true); // Start with loading true
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

  const fetchInitialResults = useCallback(async () => {
    setIsLoading(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/initial-results`);
      if (!response.ok) {
        throw new Error("Failed to fetch initial results");
      }
      const data = await response.json();
      setResults(data);
    } catch (err) {
      console.error("Error fetching initial results:", err);
      // Don't set error for initial load failures, just leave results empty
    } finally {
      setIsLoading(false);
    }
  }, []);

  useEffect(() => {
    // Load both filters and initial results when component mounts
    const loadInitialData = async () => {
      await fetchFilters();
      await fetchInitialResults();
    };
    
    loadInitialData();
  }, [fetchFilters, fetchInitialResults]);

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
        // Transform raw API results into treatment groups
        const formatted = transformApiResponseToTreatmentGroups(data.results || data);
        setResults(formatted);
        
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