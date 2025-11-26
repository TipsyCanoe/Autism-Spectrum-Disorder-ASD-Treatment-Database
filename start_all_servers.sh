#!/bin/bash
# ASD Database - Production Server Management Script
# This script manages the systemd services for production deployment

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Default environment
ENVIRONMENT=${ENVIRONMENT:-local}

echo "Managing ASD Database services for $ENVIRONMENT environment..."

# Load environment configuration for frontend setup
source "$SCRIPT_DIR/load-env.sh"

# Setup frontend environment variables and build if needed
"$SCRIPT_DIR/setup-frontend-env.sh"

if [ "$ENVIRONMENT" = "production" ]; then
    echo "Production mode: Managing systemd services..."
    
    # Build frontend for production
    echo "Building frontend..."
    cd "$SCRIPT_DIR/frontend/testing-website"
    npm run build
    
    # Restart systemd services
    echo "Restarting backend services..."
    sudo systemctl restart asd-backend.service
    
    # Reload nginx to serve new frontend build
    echo "Reloading nginx..."
    sudo nginx -t && sudo systemctl reload nginx
    
    # Check service status
    echo "Service status:"
    sudo systemctl status asd-backend.service asd-node-backend.service nginx.service --no-pager
    
else
    echo "Development mode: Starting services manually..."
    
    # Original development startup logic
    LOGS_DIR="$SCRIPT_DIR/logs"
    mkdir -p "$LOGS_DIR"
    
    # Choose Python command
    PYTHON_CMD="python3"
    command -v python3 &> /dev/null || PYTHON_CMD="python"
    
    # Memory optimization for Flask backend - reduce model cache
    export TRANSFORMERS_CACHE_MAX_ENTRIES=1
    export SENTENCE_TRANSFORMERS_HOME="$SCRIPT_DIR/.cache/sentence_transformers"
    
    # Start Frontend (with reduced memory allocation)
    echo "Starting Frontend server (port $FRONTEND_PORT)..."
    (cd "$SCRIPT_DIR/frontend/testing-website" && 
     NODE_OPTIONS="--max-old-space-size=256" PORT=$FRONTEND_PORT npm start > "$LOGS_DIR/frontend.log" 2>&1) &
    FRONTEND_PID=$!
    
    # Start Flask Backend with memory limits
    echo "Starting Backend Flask API (port $PYTHON_BACKEND_PORT)..."
    (cd "$SCRIPT_DIR/backend" && 
     source "$SCRIPT_DIR/venv/bin/activate" &&
     PYTHONUNBUFFERED=1 $PYTHON_CMD -B app.py > "$LOGS_DIR/backend_flask.log" 2>&1) &
    BACKEND_QUERY_PID=$!
    
    # Save PIDs
    echo "$FRONTEND_PID $BACKEND_QUERY_PID" > "$SCRIPT_DIR/.server_pids"
    
    echo "All servers started. PIDs: $FRONTEND_PID $BACKEND_QUERY_PID"
    echo "View logs with: tail -f $LOGS_DIR/*.log"
    echo "Stop servers with: ./stop_all_servers.sh"
    
    # Add cleanup trap
    trap 'echo "Stopping servers..."; kill $FRONTEND_PID $BACKEND_QUERY_PID 2>/dev/null; exit' INT
    
    # Keep script running
    wait $FRONTEND_PID $BACKEND_QUERY_PID
fi