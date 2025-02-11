import axios from "axios";

const search = async (selectedOptions) => {
  const response = await axios.post("/api/search", { selectedOptions });
  return response.data;
};

export default {
  search,
};
