#!/bin/bash
# Build all FIDASIM services with correct user ownership

set -e

echo "=================================="
echo "Building FIDASIM Platform"
echo "=================================="
echo ""

# Get current user's UID and GID
export USER_ID=$(id -u)
export GROUP_ID=$(id -g)

echo "Configuration:"
echo "  User ID:  $USER_ID"
echo "  Group ID: $GROUP_ID"
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

echo "=================================="
echo "Building Core Service..."
echo "=================================="
docker-compose -f docker-compose.prod.yaml build fidasim-core

echo ""
echo "=================================="
echo "Building Web Services..."
echo "=================================="
./scripts/build-web.sh

echo ""
echo "=================================="
echo "✓ Build Complete!"
echo "=================================="
echo ""
echo "Files will be owned by UID:GID = $USER_ID:$GROUP_ID"
echo ""
echo "To start services:"
echo "  USER_ID=$USER_ID GROUP_ID=$GROUP_ID docker-compose -f docker-compose.prod.yaml up -d"
echo ""
