# Nginx Configuration Notes

## Current Setup (Production)

The production nginx configuration at `/etc/nginx/sites-available/asd-db.conf` proxies:

- `/api/` → `http://127.0.0.1:5000` (Flask Python backend)
- `/jobs/` → `http://127.0.0.1:5001/` (Node.js job backend)

**Important:** The trailing slash on the Node proxy is critical!

```nginx
# Flask API - no trailing slash needed
location /api/ {
  proxy_pass http://127.0.0.1:5000;
}

# Node job backend - MUST have trailing slash to strip /jobs/ prefix
location /jobs/ {
  proxy_pass http://127.0.0.1:5001/;
}
```

## Why the Trailing Slash Matters

- **Without trailing slash:** nginx forwards the full path (e.g., `/jobs/api/run-job` → Node)
- **With trailing slash:** nginx strips the prefix (e.g., `/jobs/api/run-job` → `/api/run-job` → Node)

The Node server routes are defined as `/api/run-job`, so nginx must strip the `/jobs/` prefix.

## Frontend Configuration

The frontend uses these environment variables:

```env
REACT_APP_PYTHON_API_URL=https://star.cs.wwu.edu
REACT_APP_NODE_API_URL=https://star.cs.wwu.edu/jobs
```

The frontend code then appends:
- `/api/filters`, `/api/search`, etc. for Python API
- `/api/run-job` for Node API

## Testing

Test the endpoints:

```bash
# Python API
curl https://star.cs.wwu.edu/api/filters

# Node API
curl -X POST https://star.cs.wwu.edu/jobs/api/run-job
```

Both should return JSON responses (200 OK).
