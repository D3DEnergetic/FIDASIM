#!/bin/bash
# Deployment script for FIDASIM Web services

set -e

echo "=================================="
echo "Deploying FIDASIM Web Platform"
echo "=================================="
echo ""

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

# Check if .env exists
if [ ! -f .env ]; then
    echo "ERROR: .env file not found!"
    echo "Please create .env from .env.prod.example"
    echo ""
    echo "  cp .env.prod.example .env"
    echo "  # Edit .env with your configuration"
    echo ""
    exit 1
fi

# Load environment
echo "Loading environment from .env..."
export $(cat .env | grep -v '^#' | xargs)

# Validate critical environment variables
if [ "$SECRET_KEY" == "your-super-secret-key-change-this-in-production" ] || \
   [ "$SECRET_KEY" == "your-super-secret-key-change-this-in-production-use-openssl-rand-hex-32" ]; then
    echo "WARNING: SECRET_KEY is still set to default value!"
    echo "Generate a secure key with: openssl rand -hex 32"
    echo ""
    read -p "Continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo ""
echo "Configuration validated"
echo ""

# Set USER_ID and GROUP_ID for proper file ownership
export USER_ID=$(id -u)
export GROUP_ID=$(id -g)

echo "File ownership will use UID:GID = $USER_ID:$GROUP_ID"
echo ""

# Build services
echo "=================================="
echo "Building services..."
echo "=================================="
bash "$SCRIPT_DIR/build-all.sh"

echo ""
echo "=================================="
echo "Starting services..."
echo "=================================="

# Stop existing services
echo "Stopping existing web services..."
docker-compose -f docker-compose.prod.yaml stop fidasim-web-backend fidasim-web-frontend 2>/dev/null || true
docker-compose -f docker-compose.prod.yaml rm -f fidasim-web-backend fidasim-web-frontend 2>/dev/null || true

# Start all services
echo "Starting FIDASIM platform..."
docker-compose -f docker-compose.prod.yaml up -d

echo ""
echo "=================================="
echo "Waiting for services to be healthy..."
echo "=================================="

# Wait for backend health check
echo -n "Waiting for backend..."
max_attempts=30
attempt=0
while [ $attempt -lt $max_attempts ]; do
    if docker exec fidasim-web-backend python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/api/health').read()" 2>/dev/null; then
        echo " ✓"
        break
    fi
    echo -n "."
    sleep 2
    attempt=$((attempt + 1))
done

if [ $attempt -eq $max_attempts ]; then
    echo " ✗ Timeout"
    echo "Backend failed to start. Check logs:"
    echo "  docker logs fidasim-web-backend"
    exit 1
fi

# Wait for frontend health check
echo -n "Waiting for frontend..."
attempt=0
while [ $attempt -lt $max_attempts ]; do
    if docker exec fidasim-web-frontend wget --no-verbose --tries=1 --spider http://localhost/ 2>/dev/null; then
        echo " ✓"
        break
    fi
    echo -n "."
    sleep 2
    attempt=$((attempt + 1))
done

if [ $attempt -eq $max_attempts ]; then
    echo " ✗ Timeout"
    echo "Frontend failed to start. Check logs:"
    echo "  docker logs fidasim-web-frontend"
    exit 1
fi

echo ""
echo "=================================="
echo "✓ Deployment Complete!"
echo "=================================="
echo ""
echo "Services Status:"
docker-compose -f docker-compose.prod.yaml ps
echo ""
echo "Access Points:"
echo "  • Frontend:        http://localhost"
echo "  • Backend API:     http://localhost:8000"
echo "  • API Docs:        http://localhost:8000/api/docs"
echo "  • Core API:        http://localhost:8001"
echo "  • Preprocessor:    http://localhost:8002"
echo ""
echo "Useful Commands:"
echo "  • View logs:       docker-compose -f docker-compose.prod.yaml logs -f"
echo "  • Stop services:   docker-compose -f docker-compose.prod.yaml down"
echo "  • Restart service: docker-compose -f docker-compose.prod.yaml restart <service>"
echo ""
