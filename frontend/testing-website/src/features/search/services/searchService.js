import axios from "axios";

// Define the base URL for the backend API
const API_BASE_URL = "http://localhost:5000"; // Use the port your backend is running on

const search = async (selectedOptions, searchQuery) => {
  // Use the full URL by combining the base URL and the endpoint path
  const response = await axios.post(`${API_BASE_URL}/api/search`, {
    selectedOptions,
    searchQuery,
  });
  return response.data;
};

export default {
  search,
};
