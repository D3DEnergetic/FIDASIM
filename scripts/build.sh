#!/bin/bash

# FIDASIM Docker Build Script
# Builds all Docker images for the FIDASIM containerized application

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$PROJECT_ROOT/build_$BUILD_TIMESTAMP.log"

# Default values
BUILD_MODE="production"
SKIP_CACHE=false
PARALLEL_BUILD=false
PUSH_TO_REGISTRY=false
REGISTRY_URL=""
TAG="latest"

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
}

# Function to check prerequisites
check_prerequisites() {
    print_info "Checking prerequisites..."

    # Check Docker
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed"
        exit 1
    fi

    # Check Docker Compose
    if ! command -v docker-compose &> /dev/null; then
        print_warning "docker-compose is not installed, using docker compose"
        COMPOSE_CMD="docker compose"
    else
        COMPOSE_CMD="docker-compose"
    fi

    # Check if Docker daemon is running
    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running"
        exit 1
    fi

    # Create data directory if it doesn't exist
    if [ ! -d "$PROJECT_ROOT/data" ]; then
        print_info "Creating data directory..."
        mkdir -p "$PROJECT_ROOT/data"
    fi

    print_info "Prerequisites check completed"
}

# Function to build atomic tables if needed
prepare_atomic_tables() {
    print_info "Checking atomic tables..."

    if [ ! -f "$PROJECT_ROOT/tables/atomic_tables.h5" ]; then
        print_warning "Atomic tables not found. They will be built during Docker image creation."
        print_warning "This may take 30+ minutes on first build."

        read -p "Do you want to build atomic tables locally first? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_info "Building atomic tables locally..."
            cd "$PROJECT_ROOT"
            make tables
            cd -
        fi
    else
        print_info "Atomic tables found at tables/atomic_tables.h5"
    fi
}

# Function to build Docker images
build_images() {
    print_info "Building Docker images..."

    cd "$PROJECT_ROOT"

    # Build arguments
    BUILD_ARGS=""
    if [ "$SKIP_CACHE" = true ]; then
        BUILD_ARGS="$BUILD_ARGS --no-cache"
    fi

    # Build fidasim-core
    print_info "Building fidasim-core image..."
    if docker build $BUILD_ARGS \
        -t fidasim-core:$TAG \
        -f docker/core/Dockerfile \
        --build-arg BUILD_MODE=$BUILD_MODE \
        --progress=plain \
        . 2>&1 | tee -a "$LOG_FILE"; then
        print_info "fidasim-core image built successfully"
    else
        print_error "Failed to build fidasim-core image"
        print_error "Check the log file for details: $LOG_FILE"
        exit 1
    fi

    # Build fidasim-preprocessor
    print_info "Building fidasim-preprocessor image..."
    if docker build $BUILD_ARGS \
        -t fidasim-preprocessor:$TAG \
        -f docker/preprocessor/Dockerfile \
        --progress=plain \
        . 2>&1 | tee -a "$LOG_FILE"; then
        print_info "fidasim-preprocessor image built successfully"
    else
        print_error "Failed to build fidasim-preprocessor image"
        print_error "Check the log file for details: $LOG_FILE"
        exit 1
    fi

    # TODO: Build fidasim-web when implemented
    # print_info "Building fidasim-web image..."
    # docker build $BUILD_ARGS \
    #     -t fidasim-web:$TAG \
    #     -f docker/web/Dockerfile \
    #     . 2>&1 | tee -a "$LOG_FILE"
}

# Function to tag images for registry
tag_images() {
    if [ -n "$REGISTRY_URL" ]; then
        print_info "Tagging images for registry $REGISTRY_URL..."

        docker tag fidasim-core:$TAG $REGISTRY_URL/fidasim-core:$TAG
        docker tag fidasim-preprocessor:$TAG $REGISTRY_URL/fidasim-preprocessor:$TAG
        # docker tag fidasim-web:$TAG $REGISTRY_URL/fidasim-web:$TAG
    fi
}

# Function to push images to registry
push_images() {
    if [ "$PUSH_TO_REGISTRY" = true ] && [ -n "$REGISTRY_URL" ]; then
        print_info "Pushing images to registry $REGISTRY_URL..."

        docker push $REGISTRY_URL/fidasim-core:$TAG
        docker push $REGISTRY_URL/fidasim-preprocessor:$TAG
        # docker push $REGISTRY_URL/fidasim-web:$TAG

        print_info "Images pushed successfully"
    fi
}

# Function to verify build
verify_build() {
    print_info "Verifying build..."

    # Check if images exist
    if docker image inspect fidasim-core:$TAG &> /dev/null; then
        print_info "✓ fidasim-core:$TAG exists"
    else
        print_error "✗ fidasim-core:$TAG not found"
        exit 1
    fi

    if docker image inspect fidasim-preprocessor:$TAG &> /dev/null; then
        print_info "✓ fidasim-preprocessor:$TAG exists"
    else
        print_error "✗ fidasim-preprocessor:$TAG not found"
        exit 1
    fi

    # Show image sizes
    print_info "Image sizes:"
    docker images --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}" | grep fidasim | tee -a "$LOG_FILE"
}

# Function to run basic tests
run_tests() {
    print_info "Running basic tests..."

    # Test fidasim-core health endpoint
    print_info "Testing fidasim-core container..."
    docker run --rm fidasim-core:$TAG /app/entrypoint.sh test

    # Test fidasim-preprocessor health
    print_info "Testing fidasim-preprocessor container..."
    docker run --rm fidasim-preprocessor:$TAG python -c "import fidasim; print('FIDASIM module loaded successfully')"

    print_info "Basic tests completed"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dev|--development)
            BUILD_MODE="development"
            shift
            ;;
        --prod|--production)
            BUILD_MODE="production"
            shift
            ;;
        --no-cache)
            SKIP_CACHE=true
            shift
            ;;
        --parallel)
            PARALLEL_BUILD=true
            shift
            ;;
        --push)
            PUSH_TO_REGISTRY=true
            shift
            ;;
        --registry)
            REGISTRY_URL="$2"
            shift 2
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --test)
            RUN_TESTS=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dev, --development    Build in development mode"
            echo "  --prod, --production    Build in production mode (default)"
            echo "  --no-cache             Build without using cache"
            echo "  --parallel             Build images in parallel"
            echo "  --push                 Push images to registry"
            echo "  --registry URL         Registry URL for pushing images"
            echo "  --tag TAG              Image tag (default: latest)"
            echo "  --test                 Run tests after build"
            echo "  --help                 Show this help message"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Main execution
print_info "========================================="
print_info "FIDASIM Docker Build Script"
print_info "========================================="
print_info "Build mode: $BUILD_MODE"
print_info "Tag: $TAG"
print_info "Log file: $LOG_FILE"
print_info "========================================="

# Run build steps
check_prerequisites
prepare_atomic_tables
build_images
tag_images
push_images
verify_build

if [ "$RUN_TESTS" = true ]; then
    run_tests
fi

print_info "========================================="
print_info "Build completed successfully!"
print_info "========================================="
print_info ""
print_info "To start the services, run:"
print_info "  cd $PROJECT_ROOT"
print_info "  docker-compose up -d"
print_info ""
print_info "To view logs:"
print_info "  docker-compose logs -f"