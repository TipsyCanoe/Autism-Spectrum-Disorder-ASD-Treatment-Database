#!/bin/bash

# Get script directory and PID file
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PID_FILE="$SCRIPT_DIR/.server_pids"

echo "Stopping servers..."

# Load environment so we can stop the right ports, even if PID file is missing/stale
ENVIRONMENT=${ENVIRONMENT:-local}
if [ -f "$SCRIPT_DIR/load-env.sh" ]; then
    # shellcheck disable=SC1090
    source "$SCRIPT_DIR/load-env.sh" >/dev/null 2>&1 || true
fi

kill_by_port() {
    local port="$1"
    [ -z "$port" ] && return 0
    if command -v lsof >/dev/null 2>&1; then
        local pids
        pids=$(lsof -t -iTCP:"$port" -sTCP:LISTEN 2>/dev/null || true)
        [ -n "$pids" ] && kill $pids 2>/dev/null || true
    elif command -v fuser >/dev/null 2>&1; then
        fuser -k -n tcp "$port" >/dev/null 2>&1 || true
    fi
}

# Kill servers by PID if available
if [ -f "$PID_FILE" ]; then
    PIDS=$(cat "$PID_FILE")
    echo "Found PIDs: $PIDS"
    kill $PIDS 2>/dev/null
    sleep 1
    # Force kill if any are still running
    for PID in $PIDS; do
        ps -p "$PID" >/dev/null 2>&1 && kill -9 "$PID" 2>/dev/null
    done
    rm -f "$PID_FILE"
    # Also clean up any listeners on our configured ports
    kill_by_port "$FRONTEND_PORT"
    kill_by_port "$PYTHON_BACKEND_PORT"
    kill_by_port "$NODE_BACKEND_PORT"
    echo "Servers stopped."
else
    echo "No PID file found. Checking ports..."
    # Kill by configured ports as fallback
    kill_by_port "$FRONTEND_PORT"
    kill_by_port "$PYTHON_BACKEND_PORT"
    kill_by_port "$NODE_BACKEND_PORT"
    echo "Done."
fi