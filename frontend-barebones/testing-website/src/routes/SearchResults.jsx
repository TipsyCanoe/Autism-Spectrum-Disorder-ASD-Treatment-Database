// page for search results

import React from "react";
import "../index.css";

const SearchResults = ({ query }) => {
  return (
    <div>
      {/* Display search results based on the query */}
      <p>Results for: {query}</p>
    </div>
  );
};

export default SearchResults;