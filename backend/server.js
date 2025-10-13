require("dotenv").config(); // Load environment variables
const express = require("express");
const cors = require("cors");
const { runAPIJob } = require("./scheduler"); // Import runAPIJob

const app = express();
// Use NODE_BACKEND_PORT from environment, fallback to PORT, then 5001
const PORT = process.env.NODE_BACKEND_PORT || process.env.PORT || 5001;

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

// New Endpoint to trigger the API Job
app.post("/api/run-job", async (req, res) => {
  console.log("Received request to run API job manually...");
  try {
    const result = await runAPIJob(); // runAPIJob now returns a Promise
    console.log("API Job completed through manual trigger:", result.message);
    res.status(200).json({ 
      message: "API Job triggered successfully.", 
      details: result.message,
      output: result.output,
      stderr: result.stderr
    });
  } catch (error) {
    console.error("Error triggering API job manually:", error.message);
    res.status(500).json({ 
      message: "Failed to trigger API Job.", 
      error: error.message,
      details: error.details // Include any additional details from the error
    });
  }
});

// --- Start Server ---
app.listen(PORT, () => {
  console.log(`Backend server listening on port ${PORT}`);
});
