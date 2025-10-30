#!/bin/bash

# Debug script for troubleshooting Docker build issues

set -e

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

# Check prerequisites
print_info "Checking prerequisites..."

# Check if we're in the right directory
if [ ! -f "makefile" ] || [ ! -d "src" ]; then
    print_error "This script must be run from the FIDASIM root directory"
    exit 1
fi

# Check Docker
if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed"
    exit 1
fi

if ! docker info &> /dev/null; then
    print_error "Docker daemon is not running"
    exit 1
fi

print_info "Docker version:"
docker --version

# Check required files
print_info "Checking required files..."

REQUIRED_FILES=(
    "makefile"
    "VERSION"
    "src/fidasim.f90"
    "src/hdf5_utils.f90"
    "src/utilities.f90"
    "src/eigensystem.f90"
    "src/mpi_utils.f90"
    "deps/hdf5-1.8.16.tar.gz"
    "tables/makefile"
    "docker/core/Dockerfile"
    "docker/core/api.py"
    "docker/core/entrypoint.sh"
    "docker/preprocessor/Dockerfile"
    "docker/preprocessor/api.py"
    "docker/preprocessor/requirements.txt"
)

MISSING_FILES=()
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        MISSING_FILES+=("$file")
        print_error "Missing: $file"
    else
        print_info "✓ Found: $file"
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    print_error "Missing ${#MISSING_FILES[@]} required files"
    exit 1
fi

# Check for pre-built components
print_info ""
print_info "Checking for pre-built components..."

if [ -f "fidasim" ] && [ -x "fidasim" ]; then
    print_info "✓ Pre-built FIDASIM executable found (will be used to speed up build)"
    ls -lh fidasim
else
    print_warning "No pre-built FIDASIM executable found (will build from source)"
fi

if [ -f "tables/atomic_tables.h5" ]; then
    print_info "✓ Pre-built atomic tables found (will be used to speed up build)"
    ls -lh tables/atomic_tables.h5
else
    print_warning "No pre-built atomic tables found (will build from source - this takes 30+ minutes)"
fi

if [ -d "deps/hdf5" ] && [ -f "deps/hdf5/bin/h5fc" ]; then
    print_info "✓ Pre-built HDF5 found (will be used to speed up build)"
else
    print_warning "No pre-built HDF5 found (will build from source)"
fi

# Check data directory
print_info ""
print_info "Checking data directory..."

if [ ! -d "data" ]; then
    print_warning "Data directory not found, creating..."
    mkdir -p data
fi

print_info "Data directory: $(pwd)/data"

# Try a minimal build test
print_info ""
print_info "Testing minimal Docker build..."

# Create a test Dockerfile
cat > /tmp/test-dockerfile <<EOF
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y gfortran gcc g++ make
RUN echo "Build test successful"
EOF

if docker build -f /tmp/test-dockerfile -t test-build:latest /tmp > /dev/null 2>&1; then
    print_info "✓ Docker can build basic Ubuntu image"
    docker rmi test-build:latest > /dev/null 2>&1
else
    print_error "Docker cannot build basic Ubuntu image"
    exit 1
fi

rm -f /tmp/test-dockerfile

# Check disk space
print_info ""
print_info "Checking disk space..."
AVAILABLE_SPACE=$(df -h . | awk 'NR==2 {print $4}')
print_info "Available disk space: $AVAILABLE_SPACE"

# Show Docker system info
print_info ""
print_info "Docker system information:"
docker system df

# Try to build with more verbose output
print_info ""
print_info "Attempting to build fidasim-core with verbose output..."
print_info "This will show detailed build steps..."
echo ""

# Build with plain progress and no BuildKit
DOCKER_BUILDKIT=0 docker build \
    --progress=plain \
    --no-cache \
    -t fidasim-core:debug \
    -f docker/core/Dockerfile \
    . 2>&1 | tee build-debug.log

if [ $? -eq 0 ]; then
    print_info ""
    print_info "✓ Build successful!"
    print_info "Image created: fidasim-core:debug"
    docker images | grep fidasim-core
else
    print_error ""
    print_error "Build failed. Check build-debug.log for details"
    print_error "Common issues:"
    print_error "  1. Missing files in build context"
    print_error "  2. Network issues downloading packages"
    print_error "  3. Insufficient disk space"
    print_error "  4. Docker daemon issues"

    # Show last 20 lines of error
    print_error ""
    print_error "Last 20 lines of build log:"
    tail -20 build-debug.log

    exit 1
fi