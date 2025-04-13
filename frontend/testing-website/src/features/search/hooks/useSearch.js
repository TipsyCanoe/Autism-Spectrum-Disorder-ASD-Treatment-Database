// filepath: /home/coleoliva/senior-proj/frontend/testing-website/src/features/search/hooks/useSearch.js
import { useState, useCallback } from "react"; // Import useCallback
import searchService from "../services/searchService";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]);
  const [isLoading, setIsLoading] = useState(false); // Optional: add loading state
  const [error, setError] = useState(null); // Optional: add error state

  // Update fetchResults to accept searchQuery
  const fetchResults = useCallback(
    async (searchQuery = "") => {
      // Add searchQuery parameter
      setIsLoading(true);
      setError(null);
      try {
        // Pass both selectedOptions and searchQuery to the service
        const data = await searchService.search(selectedOptions, searchQuery);
        setResults(data);
      } catch (err) {
        console.error("Search failed:", err);
        setError(err);
        setResults([]); // Clear results on error
      } finally {
        setIsLoading(false);
      }
    },
    [selectedOptions]
  ); // Dependency array includes selectedOptions

  return {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
    isLoading, // Expose loading state
    error, // Expose error state
  };
};

export default useSearch;
