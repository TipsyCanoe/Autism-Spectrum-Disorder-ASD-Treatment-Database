import { useState } from "react";
import searchService from "../services/searchService";

const useSearch = () => {
  const [query, setQuery] = useState([]);
  const [results, setResults] = useState([]);

  const fetchResults = async (queryVector) => {
    const data = await searchService.searchVectors(queryVector);
    setResults(data);
  };

  return {
    query,
    setQuery,
    results,
    fetchResults,
  };
};

export default useSearch;
