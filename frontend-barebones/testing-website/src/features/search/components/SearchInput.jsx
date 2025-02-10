import React from "react";
import { useNavigate } from "react-router-dom";

const SearchInput = ({ query, setQuery }) => {
  const navigate = useNavigate();

  const handleSearch = () => {
    navigate("results"); // Navigate to the search results page
  };

  return (
    <div>
      <input
        type="text"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
      />
      <button onClick={handleSearch}>Search</button>
    </div>
  );
};

export default SearchInput;
