#!/bin/bash

PORTS_TO_KILL=(3000 5000 5001)
PIDS_KILLED=0

echo "Attempting to stop servers on ports: ${PORTS_TO_KILL[*]}..."
echo ""

for PORT in "${PORTS_TO_KILL[@]}"; do
    echo "Looking for process on port $PORT..."
    # Find the PID of the process using the port
    # -t: Terse output (only PID)
    # -i:$PORT: Network interface identifier (TCP/UDP on specific port)
    PID=$(lsof -t -i:$PORT)

    if [ -n "$PID" ]; then
        echo "Found process with PID $PID on port $PORT. Attempting to kill..."
        # Kill the process.
        # You can use 'kill -9 $PID' for a more forceful kill if 'kill $PID' doesn't work.
        if kill $PID; then
            echo "Successfully sent kill signal to PID $PID (port $PORT)."
            PIDS_KILLED=$((PIDS_KILLED + 1))
        else
            echo "Failed to send kill signal to PID $PID (port $PORT). It might require sudo or already be stopped."
        fi
    else
        echo "No process found running on port $PORT."
    fi
    echo ""
done

if [ "$PIDS_KILLED" -gt 0 ]; then
    echo "Finished attempting to kill processes."
else
    echo "No processes were found running on the specified ports."
fi

echo "Note: If 'npm start' opened a new terminal window/tab for the frontend server, you might need to close that manually."
