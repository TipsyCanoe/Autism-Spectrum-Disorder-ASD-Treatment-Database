# Security and Permissions Guide

This document outlines the security posture for the ASD Treatment Database deployment and the operational practices required to protect sensitive data and infrastructure.

## Access Control Overview

| Layer            | Responsible Team    | Notes |
| ---------------- | ------------------- | ----- |
| GitHub Repository| Development team    | Protected `main` branch, reviews required |
| Server (star.cs) | WWU CS Support      | SSH access limited to project maintainers |
| Database (Neon)  | Project maintainers | Role-based access via Neon console |

## Server Accounts and File Permissions

- **Deployment directory:** `/opt/asd-db`
  - Owned by application service account (e.g., `maint:maint`).
  - Only deploy user and root should have write access.
  - Use `chmod 750` on directories that should not be world-readable.
- **Systemd service files:** `/etc/systemd/system/asd-backend.service` and `/etc/systemd/system/asd-node-backend.service`
  - Owned by `root:root`, permissions `644`.
  - Edit via `sudo` only; reload daemon after changes (`sudo systemctl daemon-reload`).

## Environment Secrets

- Secrets are stored in `config/production.env` (not committed to Git).
- Ensure file permissions are restricted:
  ```bash
  sudo chmod 600 /opt/asd-db/config/production.env
  sudo chown maint:maint /opt/asd-db/config/production.env
  ```
- Never echo secrets into shell history. Use `nano` or `vi` from secure shell session.
- Rotate database credentials quarterly or when staff changes occur.

## Database Roles

- **`neondb_owner`** (provided by Neon): full privileges, used by backend.
- Create read-only roles for analytics if needed.
- Avoid sharing owner credentials with non-production systems.
- Enable SSL (`sslmode=require` included in connection string).

## Network Security

- nginx listens on ports 80/443; backends bound to localhost (127.0.0.1) only.
- Confirm firewall rules block external access to ports 5000/5001.
- Use HTTPS termination at nginx; obtain certificates via Let's Encrypt (future enhancement).

## Logging & PII Considerations

- Logs should not include user-provided email addresses or message content from Contact form.
- If sensitive data is logged accidentally, sanitize logs and rotate.
- Rotate journal logs periodically (`sudo journalctl --vacuum-time=30d`).

## Credential Rotation Checklist

1. Update credential in Neon dashboard.
2. Edit server environment file (`config/production.env`).
3. Restart backend service:
   ```bash
   sudo systemctl restart asd-backend.service
   ```
4. Verify connectivity: `curl -s https://star.cs.wwu.edu/api/initial-results`.
5. Update secret storage (GitHub Actions / future GitLab runner).

## Incident Response

- If a secret leaks:
  1. Immediately rotate the credential.
  2. Invalidate sessions / tokens if applicable.
  3. Audit access logs (`journalctl`, Neon query history).
  4. Document incident and mitigation steps in project log.
- Coordinate with WWU CS Support for server-level breaches.

## GitHub Security Settings

- Enable 2FA for all collaborators.
- Protect `main` branch with required reviews and status checks.
- Use deploy keys / GitHub Apps instead of personal tokens for automation.
- Keep GitHub Secrets restricted to necessary workflows.

Adhering to these practices keeps both patient-related research data and infrastructure secure.
