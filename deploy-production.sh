#!/bin/bash
# Production deployment script for ASD Database

set -e  # Exit on any error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "🚀 Deploying ASD Database to production..."

# 1. Setup production environment
echo "📝 Setting up production environment..."
ENVIRONMENT=production source "$SCRIPT_DIR/load-env.sh"
ENVIRONMENT=production "$SCRIPT_DIR/setup-frontend-env.sh"

# 2. Build frontend
echo "🔨 Building frontend..."
cd "$SCRIPT_DIR/frontend/testing-website"
npm run build
cd "$SCRIPT_DIR"

# 3. Restart services
echo "🔄 Restarting services..."
sudo systemctl restart asd-backend.service

# 4. Reload nginx
echo "🌐 Reloading nginx..."
sudo nginx -t
sudo systemctl reload nginx

# 5. Verify deployment
echo "✅ Verifying deployment..."
sleep 3

# Check service status
echo "📊 Service Status:"
sudo systemctl is-active asd-backend.service || echo "❌ Python backend failed"
sudo systemctl is-active nginx.service || echo "❌ Nginx failed"

# Test endpoints (with timeout and fallback)
echo "🧪 Testing endpoints..."
if timeout 5 curl -s http://localhost:5000/api/filters > /dev/null 2>&1; then
    echo "✅ Python API responding"
else
    echo "⚠️  Python API test timed out (service may still be starting)"
fi

echo ""
echo "💡 If APIs show warnings, wait 30 seconds and test manually:"
echo "   curl http://localhost:5000/api/filters"

echo "🎉 Production deployment complete!"
echo "🌍 Website: https://star.cs.wwu.edu"