#!/bin/bash

# Docker Diagnostic Script
# Helps identify Docker configuration and permission issues

set +e  # Don't exit on error so we can run all tests

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Check Docker installation
print_header "Docker Installation Check"

if command -v docker &> /dev/null; then
    print_info "Docker is installed"
    docker --version
else
    print_error "Docker is not installed"
    echo "Please install Docker from https://docs.docker.com/get-docker/"
    exit 1
fi

# Check Docker daemon
print_header "Docker Daemon Status"

if docker info > /dev/null 2>&1; then
    print_info "Docker daemon is running"
else
    print_error "Docker daemon is not running or not accessible"
    echo "Try: sudo systemctl start docker"
    echo "Or on Mac/Windows, ensure Docker Desktop is running"

    # Check if it's a permission issue
    if sudo docker info > /dev/null 2>&1; then
        print_warning "Docker works with sudo. This is a permission issue."
        echo "Fix with: sudo usermod -aG docker $USER"
        echo "Then log out and back in"
    fi
    exit 1
fi

# Check user permissions
print_header "User Permissions"

echo "Current user: $(whoami)"
echo "User groups: $(groups)"

if groups | grep -q docker; then
    print_info "User is in docker group"
else
    print_warning "User is NOT in docker group"
    echo "Fix with: sudo usermod -aG docker $USER"
    echo "Then log out and back in"
fi

# Check Docker system info
print_header "Docker System Information"

docker version
echo ""
docker system info 2>/dev/null | head -20

# Check disk space
print_header "Disk Space"

df -h / | grep -v Filesystem
df -h /var/lib/docker 2>/dev/null || df -h /var

# Check Docker storage driver
print_header "Docker Storage Configuration"

docker info 2>/dev/null | grep -E "Storage Driver|Docker Root Dir|Registry"

# Test basic Docker functionality
print_header "Docker Functionality Tests"

print_info "Test 1: Pull a small image"
if docker pull hello-world > /dev/null 2>&1; then
    print_info "✓ Can pull images"
else
    print_error "✗ Cannot pull images"
    echo "This might be a network or proxy issue"
fi

print_info "Test 2: Run a container"
if docker run --rm hello-world > /dev/null 2>&1; then
    print_info "✓ Can run containers"
else
    print_error "✗ Cannot run containers"
fi

print_info "Test 3: Build a minimal image"
cat > /tmp/Dockerfile.test <<EOF
FROM alpine:latest
RUN echo "test"
EOF

if docker build -t test-build:latest /tmp -f /tmp/Dockerfile.test > /dev/null 2>&1; then
    print_info "✓ Can build images"
    docker rmi test-build:latest > /dev/null 2>&1
else
    print_error "✗ Cannot build images"
    echo "Trying with more verbose output..."
    docker build -t test-build:latest /tmp -f /tmp/Dockerfile.test
fi
rm -f /tmp/Dockerfile.test

print_info "Test 4: Build with Ubuntu base"
cat > /tmp/Dockerfile.ubuntu <<EOF
FROM ubuntu:22.04
RUN apt-get update && echo "success"
EOF

if docker build -t test-ubuntu:latest /tmp -f /tmp/Dockerfile.ubuntu > /dev/null 2>&1; then
    print_info "✓ Can build Ubuntu-based images"
    docker rmi test-ubuntu:latest > /dev/null 2>&1
else
    print_error "✗ Cannot build Ubuntu-based images"
    echo "This is the specific issue affecting FIDASIM build"
    echo ""
    echo "Trying with detailed output..."
    DOCKER_BUILDKIT=0 docker build --progress=plain -t test-ubuntu:latest /tmp -f /tmp/Dockerfile.ubuntu 2>&1 | tail -30
fi
rm -f /tmp/Dockerfile.ubuntu

# Check for proxy settings
print_header "Network/Proxy Configuration"

if [ -n "$HTTP_PROXY" ] || [ -n "$http_proxy" ]; then
    print_warning "HTTP_PROXY is set: ${HTTP_PROXY:-$http_proxy}"
fi

if [ -n "$HTTPS_PROXY" ] || [ -n "$https_proxy" ]; then
    print_warning "HTTPS_PROXY is set: ${HTTPS_PROXY:-$https_proxy}"
fi

if [ -n "$NO_PROXY" ] || [ -n "$no_proxy" ]; then
    print_info "NO_PROXY is set: ${NO_PROXY:-$no_proxy}"
fi

# Check Docker daemon proxy settings
if [ -f /etc/systemd/system/docker.service.d/http-proxy.conf ]; then
    print_warning "Docker daemon has proxy configuration"
    cat /etc/systemd/system/docker.service.d/http-proxy.conf
fi

# Check DNS
print_header "DNS Configuration"

print_info "Testing DNS resolution..."
if nslookup archive.ubuntu.com > /dev/null 2>&1; then
    print_info "✓ DNS resolution works"
else
    print_error "✗ DNS resolution failed"
fi

print_info "Docker daemon DNS settings:"
docker info 2>/dev/null | grep -A2 "Registry"

# Check for common issues
print_header "Common Issues Check"

# Check for snap Docker
if which docker | grep -q snap; then
    print_warning "Docker is installed via snap"
    echo "Snap Docker can have permission issues. Consider installing via apt/yum instead."
fi

# Check for Docker Desktop vs Docker Engine
if [ -f /usr/local/bin/docker ]; then
    print_info "Docker Desktop detected"
elif [ -f /usr/bin/docker ]; then
    print_info "Docker Engine detected"
fi

# Check BuildKit
print_header "BuildKit Configuration"

if [ -n "$DOCKER_BUILDKIT" ]; then
    print_info "DOCKER_BUILDKIT is set to: $DOCKER_BUILDKIT"
else
    print_info "DOCKER_BUILDKIT is not set (using default)"
fi

# Final recommendations
print_header "Recommendations"

ISSUES=0

if ! docker run --rm hello-world > /dev/null 2>&1; then
    ((ISSUES++))
    print_error "Basic Docker functionality is not working"
    echo "1. Ensure Docker daemon is running"
    echo "2. Check user permissions (add to docker group)"
    echo "3. Check Docker installation"
fi

if ! docker build -t test:latest /tmp -f /tmp/Dockerfile.test > /dev/null 2>&1 2>/dev/null; then
    ((ISSUES++))
    print_error "Docker build is not working"
    echo "1. Check disk space"
    echo "2. Check Docker storage driver"
    echo "3. Try: docker system prune -a (WARNING: removes all unused images)"
fi

if [ -n "$HTTP_PROXY" ] || [ -n "$HTTPS_PROXY" ]; then
    print_warning "Proxy is configured"
    echo "Ensure Docker daemon has the same proxy settings"
    echo "See: https://docs.docker.com/config/daemon/systemd/#httphttps-proxy"
fi

if [ $ISSUES -eq 0 ]; then
    print_info "✓ Docker appears to be working correctly"
else
    print_error "Found $ISSUES issue(s) that need to be resolved"
fi

# Cleanup
docker rmi hello-world > /dev/null 2>&1 || true