import React, { useEffect } from "react";
import { useSearchParams } from "react-router-dom";
import useSearch from "../hooks/useSearch";

const SearchResults = () => {
  const [searchParams] = useSearchParams();
  const filters = JSON.parse(searchParams.get("filters") || "[]");
  const { results, fetchResults } = useSearch();

  useEffect(() => {
    if (filters.length > 0) {
      // Assuming fetchResults can handle an array of filters
      fetchResults(filters);
    }
  }, [filters, fetchResults]);

  return (
    <div>
      <h1>Search Results for: {filters.join(", ")}</h1>
      <ul>
        {results.map((result, index) => (
          <li key={index}>{result.join(", ")}</li>
        ))}
      </ul>
    </div>
  );
};

export default SearchResults;
