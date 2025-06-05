#!/bin/bash
# filepath: /home/coleoliva/senior-proj/start_all_servers.sh

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
LOGS_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOGS_DIR"

# Choose Python command
PYTHON_CMD="python3"
command -v python3 &> /dev/null || PYTHON_CMD="python"

echo "Starting servers using $PYTHON_CMD..."

# Memory optimization for Flask backend - reduce model cache
export TRANSFORMERS_CACHE_MAX_ENTRIES=1
export SENTENCE_TRANSFORMERS_HOME="$SCRIPT_DIR/.cache/sentence_transformers"

# Start Frontend (with reduced memory allocation)
echo "Starting Frontend server (port 3000)..."
(cd "$SCRIPT_DIR/frontend/testing-website" && 
 NODE_OPTIONS="--max-old-space-size=256" npm start > "$LOGS_DIR/frontend.log" 2>&1) &
FRONTEND_PID=$!

# Start Flask Backend with memory limits
echo "Starting Backend Flask API (port 5000)..."
(cd "$SCRIPT_DIR/backend" && 
 PYTHONUNBUFFERED=1 $PYTHON_CMD -B app.py > "$LOGS_DIR/backend_flask.log" 2>&1) &
BACKEND_QUERY_PID=$!

# Start Node Backend with reduced memory limits
echo "Starting Backend Node.js API (port 5001)..."
(cd "$SCRIPT_DIR/backend" && 
 NODE_OPTIONS="--max-old-space-size=256" node server.js > "$LOGS_DIR/backend_node.log" 2>&1) &
BACKEND_JOB_PID=$!

# Save PIDs
echo "$FRONTEND_PID $BACKEND_QUERY_PID $BACKEND_JOB_PID" > "$SCRIPT_DIR/.server_pids"

echo "All servers started. PIDs: $FRONTEND_PID $BACKEND_QUERY_PID $BACKEND_JOB_PID"
echo "View logs with: tail -f $LOGS_DIR/*.log"
echo "Stop servers with: ./stop_all_servers.sh"

# Add cleanup trap
trap 'echo "Stopping servers..."; kill $FRONTEND_PID $BACKEND_QUERY_PID $BACKEND_JOB_PID 2>/dev/null; exit' INT

# Keep script running
wait $FRONTEND_PID $BACKEND_QUERY_PID $BACKEND_JOB_PID