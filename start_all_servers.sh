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
    sudo systemctl restart asd-backend.service asd-node-backend.service
    
    # Reload nginx to serve new frontend build
    echo "Reloading nginx..."
    sudo nginx -t && sudo systemctl reload nginx
    
    # Check service status
    echo "Service status:"
    sudo systemctl status asd-backend.service asd-node-backend.service nginx.service --no-pager
    
else
    echo "Development mode: Starting services manually..."

    port_in_use() {
        local port="$1"
        if command -v lsof >/dev/null 2>&1; then
            lsof -iTCP:"$port" -sTCP:LISTEN -P -n >/dev/null 2>&1
            return $?
        elif command -v ss >/dev/null 2>&1; then
            ss -ltn "sport = :$port" 2>/dev/null | awk 'NR>1{found=1} END{exit(found?0:1)}'
            return $?
        elif command -v netstat >/dev/null 2>&1; then
            netstat -ltn 2>/dev/null | awk -v p=":$port" '$4 ~ p {found=1} END{exit(found?0:1)}'
            return $?
        fi
        return 1
    }

    describe_port_owner() {
        local port="$1"
        if command -v lsof >/dev/null 2>&1; then
            lsof -iTCP:"$port" -sTCP:LISTEN -P -n 2>/dev/null || true
        elif command -v ss >/dev/null 2>&1; then
            ss -ltnp "sport = :$port" 2>/dev/null || true
        elif command -v netstat >/dev/null 2>&1; then
            netstat -ltnp 2>/dev/null | awk -v p=":$port" '$4 ~ p {print}' || true
        else
            echo "(No lsof/ss/netstat available to identify the process.)"
        fi
    }

    ensure_service_started() {
        local name="$1"
        local pid="$2"
        local logfile="$3"

        sleep 1
        if ! kill -0 "$pid" >/dev/null 2>&1; then
            echo "Error: $name failed to start (PID $pid exited)."
            if [ -f "$logfile" ]; then
                echo "--- Last 50 lines of $logfile ---"
                tail -n 50 "$logfile" || true
            else
                echo "Log file not found: $logfile"
            fi
            return 1
        fi
        return 0
    }
    
    # Original development startup logic
    LOGS_DIR="$SCRIPT_DIR/logs"
    mkdir -p "$LOGS_DIR"
    
    # Choose Python command
    PYTHON_CMD="python3"
    command -v python3 &> /dev/null || PYTHON_CMD="python"
    
    # Memory optimization for Flask backend - reduce model cache
    export TRANSFORMERS_CACHE_MAX_ENTRIES=1
    export SENTENCE_TRANSFORMERS_HOME="$SCRIPT_DIR/.cache/sentence_transformers"

    if port_in_use "$FRONTEND_PORT"; then
        echo "Error: FRONTEND_PORT=$FRONTEND_PORT is already in use."
        describe_port_owner "$FRONTEND_PORT"
        echo "Tip: run ./stop_all_servers.sh or change config/${ENVIRONMENT}.env"
        exit 1
    fi

    if port_in_use "$PYTHON_BACKEND_PORT"; then
        echo "Error: PYTHON_BACKEND_PORT=$PYTHON_BACKEND_PORT is already in use."
        describe_port_owner "$PYTHON_BACKEND_PORT"
        echo "Tip: run ./stop_all_servers.sh or change config/${ENVIRONMENT}.env"
        exit 1
    fi
    
    # Start Frontend (with reduced memory allocation)
    echo "Starting Frontend server (port $FRONTEND_PORT)..."
    (cd "$SCRIPT_DIR/frontend/testing-website" && 
     NODE_OPTIONS="--max-old-space-size=256" PORT=$FRONTEND_PORT npm start > "$LOGS_DIR/frontend.log" 2>&1) &
    FRONTEND_PID=$!
    ensure_service_started "Frontend" "$FRONTEND_PID" "$LOGS_DIR/frontend.log" || {
        kill "$FRONTEND_PID" 2>/dev/null || true
        exit 1
    }
    
    # Start Flask Backend with memory limits
    echo "Starting Backend Flask API (port $PYTHON_BACKEND_PORT)..."
    (cd "$SCRIPT_DIR/services/api" && 
     if [ -f "$SCRIPT_DIR/venv/bin/activate" ]; then
        source "$SCRIPT_DIR/venv/bin/activate"
     else
        echo "Warning: venv not found at $SCRIPT_DIR/venv; using system Python." >&2
     fi &&
     PYTHONUNBUFFERED=1 $PYTHON_CMD -B app.py > "$LOGS_DIR/backend_flask.log" 2>&1) &
    BACKEND_QUERY_PID=$!
    ensure_service_started "Flask API" "$BACKEND_QUERY_PID" "$LOGS_DIR/backend_flask.log" || {
        kill "$FRONTEND_PID" 2>/dev/null || true
        kill "$BACKEND_QUERY_PID" 2>/dev/null || true
        exit 1
    }

    # Start Node.js Scheduler
    echo "Starting Node.js Scheduler..."
    (cd "$SCRIPT_DIR/services/scheduler" && 
     node server.js > "$LOGS_DIR/backend_node.log" 2>&1) &
    SCHEDULER_PID=$!
    ensure_service_started "Scheduler" "$SCHEDULER_PID" "$LOGS_DIR/backend_node.log" || {
        kill "$FRONTEND_PID" 2>/dev/null || true
        kill "$BACKEND_QUERY_PID" 2>/dev/null || true
        kill "$SCHEDULER_PID" 2>/dev/null || true
        exit 1
    }
    
    # Save PIDs
    echo "$FRONTEND_PID $BACKEND_QUERY_PID $SCHEDULER_PID" > "$SCRIPT_DIR/.server_pids"
    
    echo "All servers started. PIDs: $FRONTEND_PID $BACKEND_QUERY_PID $SCHEDULER_PID"
    echo "View logs with: tail -f $LOGS_DIR/*.log"
    echo "Stop servers with: ./stop_all_servers.sh"
    
    # Add cleanup trap
    trap 'echo "Stopping servers..."; kill $FRONTEND_PID $BACKEND_QUERY_PID $SCHEDULER_PID 2>/dev/null; exit' INT
    
    # Keep script running
    wait $FRONTEND_PID $BACKEND_QUERY_PID $SCHEDULER_PID
fi