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

  const fetchInitialResults = useCallback(async (limit = 200) => {
    setIsLoading(true);
    try {
      const url = limit ? `${API_BASE_URL}/api/initial-results?limit=${limit}` : `${API_BASE_URL}/api/initial-results`;
      const response = await fetch(url);
      if (!response.ok) {
        throw new Error("Failed to fetch initial results");
      }
      const data = await response.json();
      // Alphabetize the treatment groups for initial results
      const alphabetizedData = data.sort((a, b) => {
        const treatmentA = (a.treatment || "").toLowerCase();
        const treatmentB = (b.treatment || "").toLowerCase();
        return treatmentA.localeCompare(treatmentB);
      });
      setResults(alphabetizedData);
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
    async (query = "", limit = 200) => {
      setIsLoading(true);
      setError(null);
      
      try {
        const params = new URLSearchParams();
        params.append("query", query);
        
        if (limit) {
          params.append("limit", limit);
        }
        
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
        // API already returns data in the correct format: [{ treatment: "...", studies: [...] }]
        setResults(data);
        
      } catch (err) {
        console.error("Error fetching results:", err);
        setError("Failed to load results. Please try again.");
      } finally {
        setIsLoading(false);
      }
    },
    [selectedOptions]
  );

  return {
    selectedOptions,
    setSelectedOptions,
    results, // This will be an array of treatment objects: [{ treatment: "...", studies: [...] }]
    fetchResults,
    fetchInitialResults,
    isLoading,
    error,
    availableFilters,
    fetchFilters,
  };
};

export default useSearch;