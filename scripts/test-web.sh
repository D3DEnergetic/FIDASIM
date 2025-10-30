#!/bin/bash
# Test script for FIDASIM Web services

set -e

echo "=================================="
echo "Testing FIDASIM Web Services"
echo "=================================="
echo ""

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test function
test_endpoint() {
    local name=$1
    local url=$2
    local expected_code=${3:-200}

    echo -n "Testing $name... "

    response=$(curl -s -w "\n%{http_code}" "$url" 2>/dev/null || echo "000")
    http_code=$(echo "$response" | tail -n1)

    if [ "$http_code" == "$expected_code" ]; then
        echo -e "${GREEN}✓ OK${NC} (HTTP $http_code)"
        return 0
    else
        echo -e "${RED}✗ FAIL${NC} (HTTP $http_code, expected $expected_code)"
        return 1
    fi
}

# Check if services are running
echo "Checking if services are running..."
if ! docker ps | grep -q fidasim-web-backend; then
    echo -e "${RED}✗ Backend is not running${NC}"
    echo "Start services with: docker-compose -f docker-compose.prod.yaml up -d"
    exit 1
fi

if ! docker ps | grep -q fidasim-web-frontend; then
    echo -e "${RED}✗ Frontend is not running${NC}"
    echo "Start services with: docker-compose -f docker-compose.prod.yaml up -d"
    exit 1
fi

echo -e "${GREEN}✓ Services are running${NC}"
echo ""

# Test endpoints
echo "Testing API Endpoints:"
echo "=================================="

failed=0

# Backend health
test_endpoint "Backend health" "http://localhost:8000/api/health" || failed=$((failed + 1))

# Backend root
test_endpoint "Backend root" "http://localhost:8000/" || failed=$((failed + 1))

# Backend API docs
test_endpoint "Backend API docs" "http://localhost:8000/api/docs" || failed=$((failed + 1))

# Frontend
test_endpoint "Frontend" "http://localhost:3000/" || failed=$((failed + 1))

# Frontend index.html
test_endpoint "Frontend index" "http://localhost:3000/index.html" || failed=$((failed + 1))

echo ""
echo "Testing Core Services:"
echo "=================================="

# Core API
test_endpoint "Core API health" "http://localhost:8001/health" || failed=$((failed + 1))

# Preprocessor API
test_endpoint "Preprocessor health" "http://localhost:8002/health" || failed=$((failed + 1))

echo ""
echo "Testing Database:"
echo "=================================="

# Test PostgreSQL
echo -n "Testing PostgreSQL... "
if docker exec fidasim-postgres pg_isready -U fidasim > /dev/null 2>&1; then
    echo -e "${GREEN}✓ OK${NC}"
else
    echo -e "${RED}✗ FAIL${NC}"
    failed=$((failed + 1))
fi

# Test Redis
echo -n "Testing Redis... "
if docker exec fidasim-redis redis-cli ping > /dev/null 2>&1; then
    echo -e "${GREEN}✓ OK${NC}"
else
    echo -e "${RED}✗ FAIL${NC}"
    failed=$((failed + 1))
fi

echo ""
echo "=================================="
if [ $failed -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    echo ""
    echo "FIDASIM Web Platform is fully operational!"
    echo ""
    echo "Access the application at:"
    echo "  Frontend:  http://localhost:3000"
    echo "  Backend:   http://localhost:8000"
    echo "  API Docs:  http://localhost:8000/api/docs"
    exit 0
else
    echo -e "${RED}✗ $failed test(s) failed${NC}"
    echo ""
    echo "Check logs for details:"
    echo "  docker-compose -f docker-compose.prod.yaml logs"
    exit 1
fi
