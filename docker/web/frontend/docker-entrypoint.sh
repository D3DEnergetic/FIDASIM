#!/bin/sh
# Docker entrypoint script for FIDASIM Web Frontend
# Allows runtime environment variable injection

set -e

# Generate env-config.js with runtime environment variables
cat <<EOF > /usr/share/nginx/html/env-config.js
window.ENV = {
  VITE_API_URL: "${VITE_API_URL:-http://localhost:8000}",
  VITE_WS_URL: "${VITE_WS_URL:-ws://localhost:8000}",
  VITE_APP_NAME: "${VITE_APP_NAME:-FIDASIM Web}",
  VITE_APP_VERSION: "${VITE_APP_VERSION:-1.0.0}"
};
EOF

echo "Environment configuration generated"
echo "VITE_API_URL: ${VITE_API_URL:-http://localhost:8000}"
echo "VITE_WS_URL: ${VITE_WS_URL:-ws://localhost:8000}"

# Execute the main command
exec "$@"
