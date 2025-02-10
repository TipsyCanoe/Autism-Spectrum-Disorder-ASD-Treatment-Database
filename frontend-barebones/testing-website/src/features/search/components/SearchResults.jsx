import React, { useEffect } from "react";
import { useSearchParams } from "react-router-dom";
import useSearch from "../hooks/useSearch";

const SearchResults = () => {
  const [searchParams] = useSearchParams();
  const query = searchParams.get("q");
  const { results, fetchResults } = useSearch();

  useEffect(() => {
    if (query) {
      const queryVector = JSON.parse(query);
      fetchResults(queryVector);
    }
  }, [query, fetchResults]);

  return (
    <div>
      <h1>Search Results for: {query}</h1>
      <ul>
        {results.map((result, index) => (
          <li key={index}>{result.join(", ")}</li>
        ))}
      </ul>
    </div>
  );
};

export default SearchResults;
