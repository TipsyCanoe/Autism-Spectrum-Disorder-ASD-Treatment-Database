import React from "react";
import { Outlet } from "react-router-dom";
import SearchFilters from "../features/search/components/SearchFilters";
import SearchInput from "../features/search/components/SearchInput";
import useSearch from "../features/search/hooks/useSearch";
import "../index.css"

const SearchPage = () => {
  const { query, setQuery } = useSearch();

  return (
    <div>
      <SearchFilters />
      <SearchInput query={query} setQuery={setQuery} />
      <Outlet /> {/* This will render the child routes */}
    </div>
  );
};

export default SearchPage;
