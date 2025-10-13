#!/bin/bash

# Frontend Environment Setup Script
# This script copies the appropriate environment variables to the React .env file

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ENVIRONMENT=${ENVIRONMENT:-local}

echo "Setting up frontend environment for: $ENVIRONMENT"

# Create .env file for React from our config
FRONTEND_ENV_FILE="$SCRIPT_DIR/frontend/testing-website/.env"
CONFIG_FILE="$SCRIPT_DIR/config/${ENVIRONMENT}.env"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Extract React-specific variables from config
echo "# Auto-generated React environment file for $ENVIRONMENT" > "$FRONTEND_ENV_FILE"
echo "# Generated at: $(date)" >> "$FRONTEND_ENV_FILE"
echo "" >> "$FRONTEND_ENV_FILE"

# Extract variables that start with REACT_APP_ or are specific port variables
grep -E "^REACT_APP_" "$CONFIG_FILE" >> "$FRONTEND_ENV_FILE" 2>/dev/null || true
grep -E "^FRONTEND_PORT=" "$CONFIG_FILE" >> "$FRONTEND_ENV_FILE" 2>/dev/null || true

# Add a blank line for readability
echo "" >> "$FRONTEND_ENV_FILE"

# Also set PORT from FRONTEND_PORT if it exists
if grep -q "^FRONTEND_PORT=" "$CONFIG_FILE"; then
    FRONTEND_PORT=$(grep "^FRONTEND_PORT=" "$CONFIG_FILE" | cut -d'=' -f2)
    echo "PORT=$FRONTEND_PORT" >> "$FRONTEND_ENV_FILE"
fi

echo "Frontend .env file updated: $FRONTEND_ENV_FILE"
echo "Contents:"
cat "$FRONTEND_ENV_FILE"