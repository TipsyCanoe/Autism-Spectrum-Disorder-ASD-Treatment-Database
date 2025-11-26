#!/bin/bash
# filepath: /home/coleoliva/senior-proj/stop_all_servers.sh

# Get script directory and PID file
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PID_FILE="$SCRIPT_DIR/.server_pids"

echo "Stopping servers..."

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
    echo "Servers stopped."
else
    echo "No PID file found. Checking ports..."
    # Kill by port as fallback
    for PORT in 3000 5000; do
        PID=$(lsof -t -i:$PORT 2>/dev/null)
        [ -n "$PID" ] && kill $PID 2>/dev/null
    done
    echo "Done."
fi