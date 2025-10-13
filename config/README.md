# Environment Configuration Helper
# This file contains shared utilities for loading environment configurations

# DO NOT EDIT THE CONFIG FILES DIRECTLY IN PRODUCTION
# Instead, use environment variables or CI/CD to override values

# Available environments:
# - local: Development on your local machine
# - staging: Pre-production testing environment  
# - production: Live server environment

# How to use:
# 1. Set ENVIRONMENT variable to choose config: export ENVIRONMENT=local
# 2. Source the appropriate config: source config/${ENVIRONMENT}.env
# 3. Start your servers with the loaded configuration

# Security Note:
# - Never commit production secrets to git
# - Use GitHub Secrets or server environment variables for sensitive data
# - The DATABASE_URL and other secrets should be overridden in production