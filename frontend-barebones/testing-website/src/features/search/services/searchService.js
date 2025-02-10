import axios from "axios";

const searchVectors = async (queryVector) => {
  const response = await axios.get(
    `http://localhost:3001/search-vectors?q=${JSON.stringify(queryVector)}`
  );
  return response.data;
};

export default {
  searchVectors,
};
