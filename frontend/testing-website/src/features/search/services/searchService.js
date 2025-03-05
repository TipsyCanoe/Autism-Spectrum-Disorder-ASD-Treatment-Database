import axios from "axios";

// Update searchService.js
const search = async (selectedOptions, searchQuery) => {
  const response = await axios.post("/api/search", {
    selectedOptions,
    searchQuery,
  });
  return response.data;
};

export default {
  search,
};
