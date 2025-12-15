# Troubleshooting Guide

This guide highlights common production issues for the ASD Treatment Database and the steps to resolve them quickly. Follow the quick checklist first, then dive into the specific scenario that matches the symptoms you are seeing.

## Quick Checklist

1. **Confirm services are running**
   ```bash
   sudo systemctl status asd-backend.service
   sudo systemctl status asd-node-backend.service
   sudo systemctl status nginx.service
   ```
2. **Inspect recent logs**
   ```bash
   sudo journalctl -u asd-backend.service --since "15 minutes ago"
   sudo journalctl -u asd-node-backend.service --since "15 minutes ago"
   sudo journalctl -u nginx.service --since "15 minutes ago"
   ```
3. **Verify database connectivity**
   ```bash
   cd /opt/asd-db/backend
   source ../venv/bin/activate
   python3 -c "import psycopg2; psycopg2.connect('postgresql://...'); print('DB OK')"
   ```
4. **Check disk space and CPU**
   ```bash
   df -h
   top -c
   ```

## Common Issues and Resolutions

### Services will not start
- **Symptom:** `systemctl status` shows `failed` or crash loop.
- **Cause:** Deploy script interrupted, missing dependencies, or syntax errors.
- **Fix:**
  1. Run `sudo journalctl -xe -u asd-backend.service` to find the stack trace.
  2. If dependency missing, rerun `pip install -r services/api/requirements.txt` inside the virtualenv.
  3. Rebuild frontend if the Node service fails due to missing build: `npm install && npm run build` in `frontend/testing-website`.

### API returns empty results
- **Symptom:** `/api/search` or `/api/initial-results` returns only a few studies.
- **Cause:** The in-memory cache stored a small dataset during testing (for example when using `?limit=2`).
- **Fix:** Restart the backend to clear cache:
  ```bash
  sudo systemctl restart asd-backend.service
  ```
- **Prevention:** Avoid hitting the API with extremely small limits in production, or extend cache logic to key by limit.

### Database connection errors
- **Symptom:** Logs contain `psycopg2.OperationalError` or timeout messages.
- **Cause:** Incorrect credentials, rotated password, Neon outage, or firewall.
- **Fix:**
  1. Confirm `DATABASE_URL` in `config/production.env` matches Neon dashboard.
  2. Test connectivity from server:
     ```bash
     psql "$DATABASE_URL" -c "SELECT 1;"
     ```
  3. If using new schema, ensure migrations ran and permissions set.

### Frontend still shows old build
- **Symptom:** Website missing latest UI changes.
- **Cause:** Browser cache or frontend build not refreshed.
- **Fix:**
  1. Force rebuild: `./deploy-production.sh`.
  2. Hard-refresh browser: `Ctrl+Shift+R` (or `Cmd+Shift+R` on macOS).
  3. Confirm `/opt/asd-db/frontend/testing-website/build` timestamps update.

### nginx 404 or double `/api/api`
- **Symptom:** Browser console shows 404 with duplicated path segments.
- **Cause:** Missing trailing slash on `/jobs/` proxy in nginx config.
- **Fix:**
  1. Ensure `/etc/nginx/sites-available/asd-db.conf` has `proxy_pass http://127.0.0.1:5001/;`
  2. Test config: `sudo nginx -t`
  3. Reload: `sudo systemctl reload nginx`

### Deploy script fails mid-run
- **Symptom:** `./deploy-production.sh` stops with error, services partially restarted.
- **Cause:** Outdated dependencies, missing environment variables, build failure.
- **Fix:**
  1. Restart deploy script after addressing error.
  2. If repeated Node install failures occur, remove `frontend/testing-website/node_modules` and rerun `npm install`.
  3. Verify environment by sourcing `config/production.env` and rerunning step manually.

### Long startup times / timeouts
- **Symptom:** Health checks fail because API takes >60s to respond after restart.
- **Cause:** MedBERT model loads at startup; first request blocked until ready.
- **Fix:**
  1. Wait 60–90 seconds after restart before testing.
  2. For CI/testing, set `DISABLE_MODEL_LOADING=1`.
  3. Consider warming cache after deployment (run `/api/filters`).

### Permission denied writing logs or builds
- **Symptom:** Deploy script cannot overwrite files.
- **Cause:** File ownership changed or script run as wrong user.
- **Fix:**
  1. Ensure deployment commands run as system user with access to `/opt/asd-db`.
  2. Restore ownership if needed: `sudo chown -R maint:maint /opt/asd-db` (replace with correct user).

### Stuck Git worktree
- **Symptom:** `git pull` fails due to local changes.
- **Fix:**
  1. Stash or commit local changes.
  2. `git pull --rebase origin main`
  3. Reapply changes from stash: `git stash pop`

## When to Escalate
- Neon database experiencing outages → Contact Neon support / check status page.
- Server resource exhaustion (RAM <100MB free) → Notify infrastructure admins.
- Security incident or credential exposure → Rotate secrets immediately and notify supervisors.

Document every outage in the project log (see `logs/` directory) for future improvements.
