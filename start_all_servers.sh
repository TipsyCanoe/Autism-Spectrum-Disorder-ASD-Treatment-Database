#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "Starting all servers from $SCRIPT_DIR..."

# Determine Python command
PYTHON_CMD="python"
if command -v python3 &> /dev/null; then
    echo "Found python3, using it."
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    echo "Found python, using it."
    PYTHON_CMD="python"
else
    echo "Error: Neither python3 nor python found in PATH. Please install Python."
    exit 1
fi
echo "Using '$PYTHON_CMD' for Python scripts."
echo ""

# Start Frontend Development Server (Port 3000)
echo "Attempting to start Frontend server (localhost:3000)..."
(cd "$SCRIPT_DIR/frontend/testing-website" && npm start) &
FRONTEND_PID=$!

# Start Backend Query API Server (Python/Flask on Port 5000)
# This assumes your backend/app.py is configured to run on port 5000.
# If it needs a specific command like 'flask run --port=5000', adjust accordingly.
echo "Attempting to start Backend Query API server (localhost:5000)..."
(cd "$SCRIPT_DIR/backend" && $PYTHON_CMD app.py) &
BACKEND_QUERY_PID=$!

# Start Backend Job API Server (Node.js on Port 5001)
echo "Attempting to start Backend Job API server (localhost:5001)..."
(cd "$SCRIPT_DIR/backend" && node server.js) &
BACKEND_JOB_PID=$!

echo "-----------------------------------------------------------------"
echo "Servers are attempting to start in the background."
echo "Monitor their output in this terminal or their respective logs if configured."
echo ""
echo "To stop the servers, you can use the following PIDs:"
echo "  Frontend (npm start): $FRONTEND_PID"
echo "  Backend Query API (python app.py): $BACKEND_QUERY_PID"
echo "  Backend Job API (node server.js): $BACKEND_JOB_PID"
echo ""
echo "Example command to stop all: kill $FRONTEND_PID $BACKEND_QUERY_PID $BACKEND_JOB_PID"
echo "(If npm start opened a new window/tab, you might need to close that manually or Ctrl+C in it)"
echo "-----------------------------------------------------------------"

# You can uncomment the 'wait' commands if you want the script to
# remain active and only exit when all background jobs are done (e.g., by Ctrl+C on this script).
# wait $FRONTEND_PID
# wait $BACKEND_QUERY_PID
# wait $BACKEND_JOB_PID
