import { useState } from "react";
import searchService from "../services/searchService";

const useSearch = () => {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [results, setResults] = useState([]);

  const fetchResults = async () => {
    const data = await searchService.search(selectedOptions);
    setResults(data);
  };

  return {
    selectedOptions,
    setSelectedOptions,
    results,
    fetchResults,
  };
};

export default useSearch;
