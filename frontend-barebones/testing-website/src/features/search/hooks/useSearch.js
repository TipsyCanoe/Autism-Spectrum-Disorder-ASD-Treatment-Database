// barebones component for hooks

import { useState } from "react";

const useSearch = () => {
  // Hook logic here
  const [query, setQuery] = useState("");
  const [filters, setFilters] = useState({});

  return { query, setQuery, filters, setFilters };
};

export default useSearch;
