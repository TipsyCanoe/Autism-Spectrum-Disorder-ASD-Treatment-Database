import { useCallback, useEffect, useState } from "react"; // Add useCallback

// Define base API URL to target port 5000 for general API calls
const API_BASE_URL = "http://localhost:5000";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [availableFilters, setAvailableFilters] = useState({
    age: [],
    symptom: [],
    gender: [],
    medication: [], // Add medication to initial state
  });

  const fetchFilters = useCallback(async () => {
    // Wrap in useCallback
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
  }, []); // Empty dependency array as API_BASE_URL is constant and setAvailableFilters is stable

  useEffect(() => {
    fetchFilters();
  }, [fetchFilters]); // Now fetchFilters is stable

  const fetchResults = useCallback(
    async (query = "") => {
      // Wrap in useCallback
      setIsLoading(true);
      setError(null);
      try {
        const params = new URLSearchParams();
        params.append("query", query);
        selectedOptions.forEach((filter) => {
          // selectedOptions is a dependency
          params.append("filters", filter);
        });
        const response = await fetch(
          `${API_BASE_URL}/api/search?${params.toString()}`
        );
        if (!response.ok) {
          throw new Error("Failed to fetch search results");
        }
        const data = await response.json();
        setResults(data.results);
      } catch (err) {
        console.error("Error fetching results:", err);
        setError("Failed to load results. Please try again.");
      } finally {
        setIsLoading(false);
      }
    },
    [selectedOptions]
  ); // Add selectedOptions as a dependency

  return {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
    isLoading,
    error,
    availableFilters,
    fetchFilters, // Also return fetchFilters if SearchPage needs to call it directly, e.g., for a refresh button
  };
};

export default useSearch;
