#!/bin/bash
set -e

# Ensure uploads directory exists and has correct permissions
if [ -d "/data/uploads" ]; then
    echo "Fixing permissions for /data/uploads..."
    # Use find to fix permissions only if we own the files
    find /data/uploads -user $(id -u) -type d -exec chmod 755 {} \; 2>/dev/null || true
    find /data/uploads -user $(id -u) -type f -exec chmod 644 {} \; 2>/dev/null || true
fi

# Create uploads directory if it doesn't exist
mkdir -p /data/uploads

# Execute the main command
exec "$@"
