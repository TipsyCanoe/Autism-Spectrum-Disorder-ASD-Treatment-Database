import React, { useState } from "react";
import { useNavigate } from "react-router-dom";
import useSearch from "../hooks/useSearch";

const SearchInput = () => {
  const { setQuery } = useSearch();
  const navigate = useNavigate();
  const [filters, setFilters] = useState([]);

  const handleSearch = () => {
    setQuery(filters);
    navigate(`/search/results?filters=${JSON.stringify(filters)}`); // Navigate to the search results page with filters
  };

  const handleFilterChange = (filter) => {
    setFilters((prevFilters) =>
      prevFilters.includes(filter)
        ? prevFilters.filter((f) => f !== filter)
        : [...prevFilters, filter]
    );
  };

  return (
    <div>
      {/* Replace this with your filter checkboxes */}
      <label>
        <input
          type="checkbox"
          value="filter1"
          onChange={() => handleFilterChange("filter1")}
        />
        Filter 1
      </label>
      <label>
        <input
          type="checkbox"
          value="filter2"
          onChange={() => handleFilterChange("filter2")}
        />
        Filter 2
      </label>
      {/* Add more filters as needed */}
      <button onClick={handleSearch}>Search</button>
    </div>
  );
};

export default SearchInput;
