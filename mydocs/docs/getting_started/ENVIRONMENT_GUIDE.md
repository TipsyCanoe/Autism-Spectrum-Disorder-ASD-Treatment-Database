# Environment Configuration Usage Guide

This project now supports multiple environments (local, staging, production) with automatic configuration management.

## Quick Start

### Local Development (Default)

```bash
# Start with local configuration (default)
./start_all_servers.sh

# Or explicitly set local environment
ENVIRONMENT=local ./start_all_servers.sh
```

### Production Deployment (Recommended)

```bash
# Use the dedicated production deployment script
./deploy-production.sh

# Or manually with environment
ENVIRONMENT=production ./start_all_servers.sh
```

### Staging Environment

```bash
# Set staging environment  
ENVIRONMENT=staging ./start_all_servers.sh
```

### Production Service Management

```bash
# Start/stop/restart individual services
sudo systemctl start asd-backend.service
sudo systemctl start asd-node-backend.service
sudo systemctl restart nginx

# Check service status
sudo systemctl status asd-backend.service asd-node-backend.service nginx
```

## Environment Files

Configuration files are located in `/config/`:

- `local.env` - Local development (ports 3000, 5000, 5001)
- `staging.env` - Staging environment (ports 3001, 6000, 6001)
- `production.env` - Production environment (ports 80, 8000, 8001)

## What Gets Configured

### Backend Services

- **Python Flask API**: Uses `PYTHON_BACKEND_PORT` and `DATABASE_URL`
- **Node.js API**: Uses `NODE_BACKEND_PORT`
- **Debug Settings**: Uses `DEBUG` and `FLASK_DEBUG` flags

### Frontend

- **React Development Server**: Uses `FRONTEND_PORT`
- **API Endpoints**: Uses `REACT_APP_PYTHON_API_URL` and `REACT_APP_NODE_API_URL`

## Manual Environment Loading

You can also load environments manually:

```bash
# Load environment variables
source load-env.sh

# Check loaded configuration
echo "Python Backend: $PYTHON_BACKEND_PORT"
echo "Node Backend: $NODE_BACKEND_PORT"  
echo "Frontend: $FRONTEND_PORT"
echo "Environment: $ENVIRONMENT"
```

## Production Deployment Workflow

### Automated Deployment

Use the dedicated production deployment script for safe, automated deployments:

```bash
./deploy-production.sh
```

This script automatically:

1. Sets up production environment variables
2. Builds the React frontend
3. Restarts systemd services (Python + Node backends)
4. Reloads nginx to serve new build
5. Verifies all services are running
6. Tests API endpoints

### Manual Production Configuration

For production deployment, the `/config/production.env` is already configured for your WWU server:

```env
# Already configured for star.cs.wwu.edu
REACT_APP_PYTHON_API_URL=https://star.cs.wwu.edu/api
REACT_APP_NODE_API_URL=https://star.cs.wwu.edu/jobs
PYTHON_BACKEND_PORT=5000  # Matches nginx proxy
NODE_BACKEND_PORT=5001    # Matches nginx proxy
DEBUG=false               # Production-ready
```

## Security Notes

- Never commit production secrets to git
- Use environment variables or CI/CD secrets for sensitive data
- The DATABASE_URL and other secrets should be overridden in production
- Consider using GitHub Secrets for automated deployments

## Testing Configuration

### Local Development Testing:

```bash
# Test local environment
ENVIRONMENT=local ./start_all_servers.sh

# In another terminal, check if services are running
curl http://localhost:5000/api/filters  # Python backend
curl http://localhost:5001/api         # Node backend
curl http://localhost:3000             # Frontend
```

### Production Testing:

```bash
# Test production services
curl http://localhost:5000/api/filters  # Python backend (via systemd)
curl http://localhost:5001/api         # Node backend (via systemd)
curl https://star.cs.wwu.edu           # Frontend (via nginx)

# Check systemd service status
sudo systemctl status asd-backend.service asd-node-backend.service
```

## Troubleshooting

### Port Conflicts

If you get port conflicts, update the port numbers in the appropriate config file.

### Environment Not Loading

Make sure to source the environment:

```bash
source load-env.sh
```

### Frontend API Calls Failing

Check that `REACT_APP_PYTHON_API_URL` and `REACT_APP_NODE_API_URL` are correctly set in your environment config file.
