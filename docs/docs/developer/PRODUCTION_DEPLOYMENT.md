# Production Deployment Guide

This guide covers deploying the ASD Treatment Database to the production server at Western Washington University.

## Architecture Overview

The production system uses:

- **Nginx**: Reverse proxy + static file server (port 80/443)
- **Gunicorn**: WSGI server for Python Flask API (port 5000)
- **Node.js**: Express server for job scheduling (port 5001)
- **Systemd**: Service management for both backends
- **Neon Database**: Cloud-hosted PostgreSQL database

## Initial Server Setup (Barebones VM)

This section details how to set up a fresh Ubuntu 20.04/22.04 LTS server from scratch.

### 1. System Prerequisites

Update the system and install essential packages:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y git curl wget build-essential python3 python3-pip python3-venv python3-dev nginx
```

### 2. Install Node.js

Install Node.js (LTS version) using NodeSource:

```bash
curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
sudo apt install -y nodejs
```

### 3. Clone the Repository

Clone the repository to `/opt/asd-db`:

```bash
sudo mkdir -p /opt/asd-db
sudo chown $USER:$USER /opt/asd-db
git clone https://github.com/TipsyCanoe/Autism-Spectrum-Disorder-ASD-Treatment-Database.git /opt/asd-db
cd /opt/asd-db
```

### 4. Backend Setup (Python)

Create a virtual environment and install dependencies:

```bash
cd /opt/asd-db/services/api
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 5. Frontend Setup (React)

Install Node.js dependencies:

```bash
cd /opt/asd-db/frontend/testing-website
npm install
```

### 6. Environment Configuration

Copy the production template and edit it:

```bash
cd /opt/asd-db/config
cp production.env.template production.env
nano production.env
```

Update the `DATABASE_URL` and other settings in `production.env`.

### 7. Systemd Service Configuration

Create the backend service file `/etc/systemd/system/asd-backend.service`:

```ini
[Unit]
Description=ASD Database Python Backend
After=network.target

[Service]
User=ubuntu
WorkingDirectory=/opt/asd-db/services/api
Environment="PATH=/opt/asd-db/services/api/venv/bin"
EnvironmentFile=/opt/asd-db/config/production.env
ExecStart=/opt/asd-db/services/api/venv/bin/gunicorn --workers 3 --bind 0.0.0.0:5000 app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

Create the node scheduler service file `/etc/systemd/system/asd-node-backend.service`:

```ini
[Unit]
Description=ASD Database Node.js Scheduler
After=network.target

[Service]
User=ubuntu
WorkingDirectory=/opt/asd-db/services/scheduler
EnvironmentFile=/opt/asd-db/config/production.env
ExecStart=/usr/bin/node server.js
Restart=always

[Install]
WantedBy=multi-user.target
```

**Note:** Replace `User=ubuntu` with your actual username if different.

Enable and start the services:

```bash
sudo systemctl daemon-reload
sudo systemctl enable asd-backend.service asd-node-backend.service
sudo systemctl start asd-backend.service asd-node-backend.service
```

### 8. Nginx Configuration

Remove the default configuration:

```bash
sudo rm /etc/nginx/sites-enabled/default
```

Create the site configuration `/etc/nginx/sites-available/asd-db.conf`:

```nginx
server {
    listen 80;
    server_name star.cs.wwu.edu;  # Replace with your domain

    root /opt/asd-db/frontend/testing-website/build;
    index index.html;

    # Serve React app
    location / {
      try_files $uri $uri/ /index.html;
    }

    # Flask API
    location /api/ {
      proxy_pass http://127.0.0.1:5000;
      proxy_set_header Host $host;
      proxy_set_header X-Real-IP $remote_addr;
    }

    # Node job backend
    location /jobs/ {
      proxy_pass http://127.0.0.1:5001/;
      proxy_http_version 1.1;
      proxy_set_header Upgrade $http_upgrade;
      proxy_set_header Connection "upgrade";
    }

    proxy_read_timeout 90s;
}
```

Enable the site and restart Nginx:

```bash
sudo ln -s /etc/nginx/sites-available/asd-db.conf /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

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

For detailed nginx configuration documentation, see [Nginx Configuration](NGINX_CONFIGURATION.md).

### Application Configuration

- `/opt/asd-db/config/production.env` - Production environment variables

## Automated Deployment (CI/CD)

The project supports automated deployment via GitHub Actions. When enabled, pushing to the `main` branch will automatically:

1. Run tests
2. Deploy to production if tests pass
3. Restart services and reload nginx

For detailed CI/CD documentation, see [CI/CD Pipeline](CI_CD_PIPELINE.md).

Note: Auto-deployment is currently disabled to isolate production during active development.

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
