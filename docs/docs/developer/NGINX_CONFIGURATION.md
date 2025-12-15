# Nginx Configuration

This document explains the nginx reverse proxy configuration for the ASD Treatment Database production deployment.

## Overview

Nginx serves multiple roles in the production architecture:

1. Static file server for the React frontend build
2. Reverse proxy for the Python Flask API
3. Reverse proxy for the Node.js job scheduler API
4. SSL/TLS termination (if configured)

## Configuration Location

The main configuration file is located at:

```
/etc/nginx/sites-available/asd-db.conf
```

It is symlinked to `/etc/nginx/sites-enabled/asd-db.conf` to enable the site.

## Current Configuration

```nginx
server {
    listen 80;
    server_name star.cs.wwu.edu;

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

    # Increase timeout for long requests
    proxy_read_timeout 90s;
}
```

## URL Routing Explained

### Frontend (React)

All requests to the root path are served from the static build directory:

- `https://star.cs.wwu.edu/` → `/opt/asd-db/frontend/testing-website/build/index.html`
- `https://star.cs.wwu.edu/search` → `index.html` (React Router handles client-side routing)

The `try_files` directive ensures that all routes fallback to `index.html`, allowing React Router to work properly.

### Python API (Flask)

API requests starting with `/api/` are proxied to the Flask backend:

- `https://star.cs.wwu.edu/api/filters` → `http://127.0.0.1:5000/api/filters`
- `https://star.cs.wwu.edu/api/search` → `http://127.0.0.1:5000/api/search`
- `https://star.cs.wwu.edu/api/initial-results` → `http://127.0.0.1:5000/api/initial-results`

**Note:** The `/api/` prefix is preserved and forwarded to Flask because `proxy_pass` does not have a trailing slash.

### Node.js API (Job Scheduler)

Job-related requests starting with `/jobs/` are proxied to the Node backend:

- `https://star.cs.wwu.edu/jobs/api/run-job` → `http://127.0.0.1:5001/api/run-job`

**Important:** The trailing slash in `proxy_pass http://127.0.0.1:5001/;` strips the `/jobs/` prefix before forwarding to the Node server. This is critical because the Node server defines routes as `/api/run-job`, not `/jobs/api/run-job`.

## Trailing Slash Behavior

The trailing slash in `proxy_pass` directives controls path rewriting:

### Without Trailing Slash

```nginx
location /api/ {
  proxy_pass http://127.0.0.1:5000;
}
```

**Behavior:** Full path is preserved
- Request: `/api/filters`
- Forwarded to: `http://127.0.0.1:5000/api/filters`

### With Trailing Slash

```nginx
location /jobs/ {
  proxy_pass http://127.0.0.1:5001/;
}
```

**Behavior:** Location prefix is stripped
- Request: `/jobs/api/run-job`
- Forwarded to: `http://127.0.0.1:5001/api/run-job` (the `/jobs/` prefix is removed)

## Frontend Environment Configuration

The React frontend is configured to use these base URLs:

```env
REACT_APP_PYTHON_API_URL=https://star.cs.wwu.edu
REACT_APP_NODE_API_URL=https://star.cs.wwu.edu/jobs
```

The frontend code then appends specific endpoints:

**Python API:**
- Appends `/api/filters` → `https://star.cs.wwu.edu/api/filters`
- Appends `/api/search` → `https://star.cs.wwu.edu/api/search`

**Node API:**
- Appends `/api/run-job` → `https://star.cs.wwu.edu/jobs/api/run-job`

## Modifying the Configuration

If you need to change the nginx configuration:

1. Edit the configuration file:
   ```bash
   sudo nano /etc/nginx/sites-available/asd-db.conf
   ```

2. Test the configuration:
   ```bash
   sudo nginx -t
   ```

3. If the test passes, reload nginx:
   ```bash
   sudo systemctl reload nginx
   ```

4. If the test fails, review error messages and fix syntax issues before reloading.

## Testing Endpoints

Verify that all proxies are working correctly:

```bash
# Test Python API
curl https://star.cs.wwu.edu/api/filters

# Test Node API
curl -X POST https://star.cs.wwu.edu/jobs/api/run-job

# Test frontend
curl https://star.cs.wwu.edu
```

All should return HTTP 200 status codes with appropriate content.

## Common Issues

### 502 Bad Gateway

**Symptoms:** Nginx returns a 502 error when accessing API endpoints.

**Causes:**
- Backend service is not running
- Backend service is listening on the wrong port
- Firewall blocking localhost connections

**Solutions:**
```bash
# Check if services are running
sudo systemctl status asd-backend.service
sudo systemctl status asd-node-backend.service

# Check if services are listening on correct ports
sudo netstat -tlnp | grep -E '5000|5001'

# Restart services
sudo systemctl restart asd-backend.service asd-node-backend.service
```

### 404 Not Found on API Routes

**Symptoms:** API requests return 404, but direct service access works.

**Causes:**
- Incorrect `proxy_pass` configuration
- Missing or incorrect trailing slash
- Route not defined in backend application

**Solutions:**
1. Verify the backend route exists:
   ```bash
   # For Python API
   grep -r "route.*filters" /opt/asd-db/services/api/app.py
   
   # For Node API
   grep -r "post.*run-job" /opt/asd-db/services/scheduler/server.js
   ```

2. Test backend directly:
   ```bash
   curl http://localhost:5000/api/filters
   curl -X POST http://localhost:5001/api/run-job
   ```

3. If direct access works but proxied access fails, check the `proxy_pass` trailing slash.

### React Router Not Working

**Symptoms:** Direct navigation to `/search` returns 404.

**Cause:** Missing `try_files` directive or incorrect configuration.

**Solution:** Ensure the location block has:
```nginx
location / {
  try_files $uri $uri/ /index.html;
}
```

This ensures all unmatched routes fallback to `index.html`, allowing React Router to handle routing.

## Performance Tuning

### Increase Timeout for Long Requests

If API requests take longer than the default timeout:

```nginx
proxy_read_timeout 120s;
proxy_connect_timeout 120s;
proxy_send_timeout 120s;
```

### Enable Gzip Compression

To reduce bandwidth usage:

```nginx
gzip on;
gzip_vary on;
gzip_min_length 1000;
gzip_types text/plain text/css application/json application/javascript text/xml application/xml;
```

### Static Asset Caching

To improve frontend performance:

```nginx
location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg|woff|woff2|ttf|eot)$ {
    expires 1y;
    add_header Cache-Control "public, immutable";
}
```

## Security Considerations

### Rate Limiting

Protect against abuse by limiting request rates:

```nginx
limit_req_zone $binary_remote_addr zone=api_limit:10m rate=10r/s;

location /api/ {
    limit_req zone=api_limit burst=20 nodelay;
    proxy_pass http://127.0.0.1:5000;
}
```

### SSL/TLS Configuration

For production deployments, configure SSL:

```nginx
server {
    listen 443 ssl http2;
    server_name star.cs.wwu.edu;

    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;
    
    # Strong SSL settings
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers HIGH:!aNULL:!MD5;
    ssl_prefer_server_ciphers on;
    
    # ... rest of configuration
}

# Redirect HTTP to HTTPS
server {
    listen 80;
    server_name star.cs.wwu.edu;
    return 301 https://$server_name$request_uri;
}
```

## Backup Configuration

Before making changes, always backup the current configuration:

```bash
sudo cp /etc/nginx/sites-available/asd-db.conf /etc/nginx/sites-available/asd-db.conf.backup
```

To restore from backup:

```bash
sudo cp /etc/nginx/sites-available/asd-db.conf.backup /etc/nginx/sites-available/asd-db.conf
sudo nginx -t
sudo systemctl reload nginx
```
