require("dotenv").config(); // Keep for PORT, maybe other future config
const express = require("express");
const cors = require("cors");
require("./scheduler");
// const { Pool } = require('pg'); // Remove or comment out DB connection for now

const app = express();
const PORT = process.env.PORT || 5001;

// --- Database Configuration ---
// Comment out or remove DB pool setup for now
// const pool = new Pool({ ... });

// --- Middleware ---
app.use(cors());
app.use(express.json());

// --- Routes ---

app.get("/api", (req, res) => {
  res.json({ message: "Backend server is running!" });
});

// Search Endpoint - Bare Minimum Placeholder
app.post("/api/search", async (req, res) => {
  // Destructure data from the frontend request body
  const { selectedOptions, searchQuery } = req.body;

  console.log("Received search request (placeholder):");
  console.log("Selected Options:", selectedOptions);
  console.log("Search Query:", searchQuery);

  // --- Placeholder Logic ---
  // Instead of querying DB or LLM, just send a success response.
  // You can optionally send back some dummy data or echo the input.
  try {
    // Simulate finding some results based on input (optional)
    const dummyResults = [
      {
        id: 1,
        title: `Placeholder result for query: "${searchQuery}"`,
        description: "This is a dummy result. LLM/DB logic is pending.",
        tags: selectedOptions || [],
        url: "#",
      },
      {
        id: 2,
        title: "Another placeholder",
        description: `Filters received: ${
          selectedOptions?.join(", ") || "None"
        }`,
        tags: ["dummy", "placeholder"],
        url: "#",
      },
    ];

    // Send back the dummy data as JSON
    res.json(dummyResults);
  } catch (error) {
    // Basic error handling in case something unexpected happens
    console.error("Error in placeholder search endpoint:", error);
    res
      .status(500)
      .json({ message: "Error in placeholder search", error: error.message });
  }
});

// --- Start Server ---
app.listen(PORT, () => {
  console.log(`Backend server listening on port ${PORT}`);
});
