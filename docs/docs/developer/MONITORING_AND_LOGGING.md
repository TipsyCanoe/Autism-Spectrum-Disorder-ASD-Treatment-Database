# Monitoring and Logging Guide

This guide explains how to observe system health, gather logs, and create actionable alerts for the ASD Treatment Database deployment.

## Key Components to Monitor

| Component         | Location / Port | Health Indicators |
| ----------------- | --------------- | ----------------- |
| nginx             | 80 / 443        | Serving static assets, proxy responses <200ms |
| Python backend    | 127.0.0.1:5000  | API uptime, response times <2s, memory usage <1.5GB |
| Node job backend  | 127.0.0.1:5001  | Job executions succeed, queue empty |
| Neon PostgreSQL   | Cloud           | Query latency <500ms, connection count <20 |

## Log Locations

- **nginx access/error logs:** `/var/log/nginx/access.log`, `/var/log/nginx/error.log`
- **Python backend:** `journalctl -u asd-backend.service`
- **Node backend:** `journalctl -u asd-node-backend.service`
- **System messages:** `journalctl -u nginx.service`, `/var/log/syslog`

### Tail Logs in Real Time
```bash
sudo journalctl -u asd-backend.service -f
sudo journalctl -u asd-node-backend.service -f
sudo tail -f /var/log/nginx/access.log
```

## Routine Health Checks

### Daily (or per deploy)
1. **Ping primary API endpoints**
   ```bash
   curl -s https://star.cs.wwu.edu/api/filters | jq '.available_filters.medication | length'
   curl -s "https://star.cs.wwu.edu/api/search?query=autism" | jq 'length'
   ```
2. **Verify frontend build timestamp**
   ```bash
   stat /opt/asd-db/frontend/testing-website/build/index.html
   ```
3. **Check system resource usage**
   ```bash
   free -h
   df -h /opt
   top -b -n1 | head -20
   ```

### Weekly
- Review `journalctl` for warnings or stack traces.
- Confirm scheduled jobs ran (if automated data updates re-enabled).
- Ensure SSL certificates (if used) have >30 days validity.

## Metrics to Track

- **Response latency:** Add simple curl timing when debugging:
  ```bash
  curl -o /dev/null -s -w 'Total: %{time_total}s\n' https://star.cs.wwu.edu/api/initial-results
  ```
- **Memory usage:** `systemctl status` shows RSS; for detailed view use `ps -o pid,ppid,cmd,%mem,%cpu -p $(pgrep -f gunicorn)`.
- **Open connections:** `ss -tulwn | grep 5000`
- **Database performance:** Use Neon dashboard query statistics; consider enabling slow query logs.

## Alerting Suggestions

Automation is minimal today, but you can layer simple alerts:

- **Simple cron-based health check:** Create a cron job that hits `/api/filters` every 10 minutes and emails on failure.
- **Log-based alerts:** Use `journalctl --since "10 minutes ago" | grep -i error` in a cronjob and notify on matches.
- **Neon alerts:** Configure alert thresholds for connection saturation and storage usage in the Neon console.

## Log Retention & Rotation

- Journald defaults to limited retention; enforce size limit:
  ```bash
  sudo journalctl --vacuum-size=500M
  ```
- nginx logs rotate via `/etc/logrotate.d/nginx`; ensure rotation is active to avoid disk pressure.
- Archive important incident logs under `/opt/asd-db/logs/` with timestamped filenames.

## Observability Roadmap (Future Enhancements)

1. **Centralized logging:** Forward journald + nginx logs to a remote syslog or ELK stack.
2. **Metrics dashboard:** Deploy Prometheus node exporter + Grafana for CPU, memory, and response time charts.
3. **Synthetic monitoring:** Use external uptime monitor (e.g., UptimeRobot) for `https://star.cs.wwu.edu`.
4. **Structured application logs:** Adopt JSON logging in Flask/Node for easier parsing.

Consistent monitoring and log review ensure issues are caught early and resolved before they impact users.
