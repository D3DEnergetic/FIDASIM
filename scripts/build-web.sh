#!/bin/bash
# Build script for FIDASIM Web services

set -e

echo "=================================="
echo "Building FIDASIM Web Services"
echo "=================================="
echo ""

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

# Source environment variables if .env exists
if [ -f .env ]; then
    echo "Loading environment from .env..."
    export $(cat .env | grep -v '^#' | xargs)
fi

# Default values
VITE_API_URL=${VITE_API_URL:-http://localhost:8000}
VITE_WS_URL=${VITE_WS_URL:-ws://localhost:8000}
SECRET_KEY=${SECRET_KEY:-your-super-secret-key-change-this-in-production}

echo "Configuration:"
echo "  API URL: $VITE_API_URL"
echo "  WS URL: $VITE_WS_URL"
echo ""

# Build Backend
echo "=================================="
echo "Building Backend..."
echo "=================================="
docker build \
    -t fidasim-web-backend:latest \
    -f docker/web/backend/Dockerfile \
    docker/web/backend

echo ""
echo "✓ Backend build complete"
echo ""

# Build Frontend
echo "=================================="
echo "Building Frontend..."
echo "=================================="
docker build \
    -t fidasim-web-frontend:latest \
    -f docker/web/frontend/Dockerfile \
    --build-arg VITE_API_URL="$VITE_API_URL" \
    --build-arg VITE_WS_URL="$VITE_WS_URL" \
    --build-arg VITE_APP_NAME="FIDASIM Web" \
    --build-arg VITE_APP_VERSION="1.0.0" \
    docker/web/frontend

echo ""
echo "✓ Frontend build complete"
echo ""

# Show built images
echo "=================================="
echo "Built Images:"
echo "=================================="
docker images | grep -E "fidasim-web|REPOSITORY"

echo ""
echo "=================================="
echo "✓ Build Complete!"
echo "=================================="
echo ""
echo "To start the services:"
echo "  docker-compose -f docker-compose.prod.yaml up -d fidasim-web-backend fidasim-web-frontend"
echo ""
