import { useState, useEffect } from "react";

// Define base API URL to target port 5000
const API_BASE_URL = "http://localhost:5000";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [availableFilters, setAvailableFilters] = useState({
    age: [],
    symptom: [],
    gender: []
  });

  useEffect(() => {
    fetchFilters();
  }, []);

  const fetchFilters = async () => {
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
  };

  const fetchResults = async (query = "") => {
    setIsLoading(true);
    setError(null);
    try {
      const params = new URLSearchParams();
      params.append("query", query);
      selectedOptions.forEach(filter => {
        params.append("filters", filter);
      });
      const response = await fetch(`${API_BASE_URL}/api/search?${params.toString()}`);
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
  };

  return {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
    isLoading,
    error,
    availableFilters
  };
};

export default useSearch;