// essentially the search page barebones

import React from "react";
import useSearch from '../features/search/hooks/useSearch';

const SearchPage = () => {
  const { query, setQuery } = useSearch();

  return (
    <div>
      <input
        type="text"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
      />
      <button>Search</button>
    </div>
  );
};

export default SearchPage;