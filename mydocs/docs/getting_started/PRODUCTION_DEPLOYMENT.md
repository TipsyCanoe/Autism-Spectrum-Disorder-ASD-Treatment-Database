# Production Deployment Guide

This guide covers deploying the ASD Treatment Database to the production server at Western Washington University.

## Architecture Overview

The production system uses:

- **Nginx**: Reverse proxy + static file server (port 80/443)
- **Gunicorn**: WSGI server for Python Flask API (port 5000)
- **Node.js**: Express server for job scheduling (port 5001)
- **Systemd**: Service management for both backends
- **Neon Database**: Cloud-hosted PostgreSQL database

## Deployment Methods

### Automated Deployment (Recommended)

```bash
./deploy-production.sh
```

This single command handles everything safely:

1. Sets production environment
2. Builds React frontend
3. Restarts backend services
4. Reloads nginx
5. Verifies deployment

### Manual Deployment

If you need more control:

```bash
# 1. Set environment and build frontend
ENVIRONMENT=production ./setup-frontend-env.sh
cd frontend/testing-website && npm run build && cd ../..

# 2. Restart services
sudo systemctl restart asd-backend.service
sudo systemctl restart asd-node-backend.service

# 3. Reload nginx
sudo nginx -t && sudo systemctl reload nginx
```

## Service Management

### Check Service Status

```bash
sudo systemctl status asd-backend.service      # Python Flask API
sudo systemctl status asd-node-backend.service # Node.js job scheduler
sudo systemctl status nginx.service            # Web server
```

### View Service Logs

```bash
sudo journalctl -fu asd-backend.service      # Python backend logs
sudo journalctl -fu asd-node-backend.service # Node backend logs
sudo journalctl -fu nginx.service            # Nginx logs
```

### Start/Stop Services

```bash
# Stop services
sudo systemctl stop asd-backend.service asd-node-backend.service

# Start services  
sudo systemctl start asd-backend.service asd-node-backend.service

# Restart services
sudo systemctl restart asd-backend.service asd-node-backend.service
```

## Configuration Files

### Systemd Services

- `/etc/systemd/system/asd-backend.service` - Python Flask service
- `/etc/systemd/system/asd-node-backend.service` - Node.js service

### Nginx Configuration

- `/etc/nginx/sites-available/asd-db.conf` - Main site config
- `/etc/nginx/sites-enabled/asd-db.conf` - Enabled site (symlink)

### Application Configuration

- `/opt/asd-db/config/production.env` - Production environment variables

## Troubleshooting

### Services Won't Start

```bash
# Check detailed error logs
sudo journalctl -xe -u asd-backend.service
sudo journalctl -xe -u asd-node-backend.service

# Check nginx configuration
sudo nginx -t

# Verify file permissions
ls -la /opt/asd-db/
```

### Database Connection Issues

```bash
# Test database connectivity
cd /opt/asd-db/backend
source ../venv/bin/activate  
python3 -c "import psycopg2; print('Database connection OK')"
```

### Frontend Not Loading

```bash
# Check if build exists
ls -la /opt/asd-db/frontend/testing-website/build/

# Rebuild if needed
cd /opt/asd-db/frontend/testing-website
npm run build

# Check nginx is serving correctly
curl -I https://star.cs.wwu.edu
```

## Health Checks

Test all components:

```bash
# API endpoints
curl -f https://star.cs.wwu.edu/api/filters   # Python API via nginx
curl -f https://star.cs.wwu.edu/jobs/api      # Node API via nginx

# Direct backend access (internal)
curl -f http://localhost:5000/api/filters     # Python API direct
curl -f http://localhost:5001/api             # Node API direct

# Frontend
curl -f https://star.cs.wwu.edu               # Main website
```

## Security Notes

- Services run as `maint` user (non-root)
- Database credentials stored in environment files
- Nginx handles SSL/TLS termination
- Backend services only accessible via localhost
- All external traffic goes through nginx reverse proxy

## Rollback Procedure

If deployment fails:

```bash
# 1. Restore previous frontend build (if you have backup)
# 2. Restart services to clear any issues
sudo systemctl restart asd-backend.service asd-node-backend.service

# 3. Check logs for issues
sudo journalctl -xe -u asd-backend.service

# 4. Reload nginx to ensure consistency
sudo systemctl reload nginx
```
