#!/bin/bash

# Environment Configuration Loader
# This script loads the appropriate environment configuration

# Default to local environment if not specified
ENVIRONMENT=${ENVIRONMENT:-local}

# Validate environment
case "$ENVIRONMENT" in
    local|staging|production)
        echo "Loading $ENVIRONMENT environment configuration..."
        ;;
    *)
        echo "Error: Invalid environment '$ENVIRONMENT'. Use 'local', 'staging', or 'production'"
        exit 1
        ;;
esac

# Load the environment configuration
CONFIG_FILE="$(dirname "$0")/config/${ENVIRONMENT}.env"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Export all variables from the config file
set -a  # Automatically export all variables
source "$CONFIG_FILE"
set +a  # Stop automatically exporting

echo "Environment configuration loaded from: $CONFIG_FILE"
echo "Python Backend Port: $PYTHON_BACKEND_PORT"
echo "Node Backend Port: $NODE_BACKEND_PORT" 
echo "Frontend Port: $FRONTEND_PORT"
echo "Environment: $ENVIRONMENT"