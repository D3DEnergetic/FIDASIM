#!/bin/bash

# FIDASIM Docker Deployment Script
# Manages deployment of FIDASIM containerized services

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DEPLOYMENT_MODE="development"
ACTION="up"

# Docker Compose command detection
if command -v docker-compose &> /dev/null; then
    COMPOSE_CMD="docker-compose"
else
    COMPOSE_CMD="docker compose"
fi

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_status() {
    echo -e "${BLUE}[STATUS]${NC} $1"
}

# Function to check environment
check_environment() {
    print_info "Checking deployment environment..."

    # Check if .env file exists
    if [ ! -f "$PROJECT_ROOT/.env" ]; then
        print_warning ".env file not found. Creating from template..."
        create_env_file
    fi

    # Check if data directory exists
    if [ ! -d "$PROJECT_ROOT/data" ]; then
        print_info "Creating data directory..."
        mkdir -p "$PROJECT_ROOT/data"
    fi

    # Check if Docker is running
    if ! docker info &> /dev/null; then
        print_error "Docker is not running"
        exit 1
    fi

    # Check if images are built
    if ! docker image inspect fidasim-core:latest &> /dev/null; then
        print_warning "fidasim-core image not found"
        read -p "Do you want to build images first? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            "$SCRIPT_DIR/build.sh"
        else
            print_error "Cannot deploy without images"
            exit 1
        fi
    fi
}

# Function to create .env file
create_env_file() {
    cat > "$PROJECT_ROOT/.env" <<EOF
# FIDASIM Docker Environment Variables

# Database Configuration
POSTGRES_DB=fidasim
POSTGRES_USER=fidasim
POSTGRES_PASSWORD=$(openssl rand -base64 32)

# Redis Configuration
REDIS_PASSWORD=$(openssl rand -base64 32)

# Application Configuration
SECRET_KEY=$(openssl rand -hex 32)
DEBUG=false
LOG_LEVEL=INFO

# Resource Limits
MAX_CORES=32
MAX_MEMORY=8G
MAX_JOBS=100

# Network Configuration
API_HOST=0.0.0.0
CORE_API_PORT=8001
PREPROCESSOR_API_PORT=8002
WEB_API_PORT=8000

# Storage
DATA_DIR=/data

# Authentication (change these!)
ADMIN_USERNAME=admin
ADMIN_PASSWORD=changeme123

# Optional: External Services
# LDAP_SERVER=
# OAUTH_PROVIDER=
# S3_BUCKET=
EOF

    print_info ".env file created. Please update passwords before production deployment!"
}

# Function to start services
start_services() {
    print_info "Starting FIDASIM services..."

    cd "$PROJECT_ROOT"

    if [ "$DEPLOYMENT_MODE" = "development" ]; then
        print_info "Starting in development mode..."
        $COMPOSE_CMD -f docker-compose.yaml -f docker-compose.dev.yaml up -d
    else
        print_info "Starting in production mode..."
        $COMPOSE_CMD -f docker-compose.yaml up -d
    fi

    # Wait for services to be healthy
    print_info "Waiting for services to be healthy..."
    sleep 5

    # Check service health
    check_service_health
}

# Function to stop services
stop_services() {
    print_info "Stopping FIDASIM services..."

    cd "$PROJECT_ROOT"

    if [ "$DEPLOYMENT_MODE" = "development" ]; then
        $COMPOSE_CMD -f docker-compose.yaml -f docker-compose.dev.yaml down
    else
        $COMPOSE_CMD -f docker-compose.yaml down
    fi

    print_info "Services stopped"
}

# Function to restart services
restart_services() {
    print_info "Restarting FIDASIM services..."
    stop_services
    sleep 2
    start_services
}

# Function to check service health
check_service_health() {
    print_info "Checking service health..."

    # Check fidasim-core
    if curl -s -f http://localhost:8001/health > /dev/null 2>&1; then
        print_status "✓ fidasim-core is healthy"
    else
        print_warning "✗ fidasim-core is not responding"
    fi

    # Check fidasim-preprocessor
    if curl -s -f http://localhost:8002/health > /dev/null 2>&1; then
        print_status "✓ fidasim-preprocessor is healthy"
    else
        print_warning "✗ fidasim-preprocessor is not responding"
    fi

    # Check PostgreSQL
    if docker exec fidasim-postgres pg_isready -U fidasim > /dev/null 2>&1; then
        print_status "✓ PostgreSQL is healthy"
    else
        print_warning "✗ PostgreSQL is not responding"
    fi

    # Check Redis
    if docker exec fidasim-redis redis-cli ping > /dev/null 2>&1; then
        print_status "✓ Redis is healthy"
    else
        print_warning "✗ Redis is not responding"
    fi

    # TODO: Check fidasim-web when implemented
}

# Function to show logs
show_logs() {
    SERVICE="$1"

    if [ -z "$SERVICE" ]; then
        print_info "Showing logs for all services..."
        cd "$PROJECT_ROOT"
        $COMPOSE_CMD logs -f
    else
        print_info "Showing logs for $SERVICE..."
        cd "$PROJECT_ROOT"
        $COMPOSE_CMD logs -f "$SERVICE"
    fi
}

# Function to show status
show_status() {
    print_info "FIDASIM Service Status"
    print_info "======================"

    cd "$PROJECT_ROOT"
    $COMPOSE_CMD ps

    echo ""
    check_service_health

    echo ""
    print_info "Resource Usage:"
    docker stats --no-stream --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}" | grep fidasim || true
}

# Function to backup database
backup_database() {
    BACKUP_DIR="$PROJECT_ROOT/backups"
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    BACKUP_FILE="$BACKUP_DIR/fidasim_db_$TIMESTAMP.sql"

    print_info "Backing up database..."

    mkdir -p "$BACKUP_DIR"

    docker exec fidasim-postgres pg_dump -U fidasim fidasim > "$BACKUP_FILE"

    if [ $? -eq 0 ]; then
        print_info "Database backed up to: $BACKUP_FILE"
        # Compress the backup
        gzip "$BACKUP_FILE"
        print_info "Backup compressed: ${BACKUP_FILE}.gz"
    else
        print_error "Database backup failed"
        exit 1
    fi
}

# Function to restore database
restore_database() {
    BACKUP_FILE="$1"

    if [ -z "$BACKUP_FILE" ]; then
        print_error "Please provide a backup file path"
        exit 1
    fi

    if [ ! -f "$BACKUP_FILE" ]; then
        print_error "Backup file not found: $BACKUP_FILE"
        exit 1
    fi

    print_warning "This will replace the current database. Are you sure? (y/n)"
    read -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Restore cancelled"
        exit 0
    fi

    print_info "Restoring database from: $BACKUP_FILE"

    # Decompress if needed
    if [[ "$BACKUP_FILE" == *.gz ]]; then
        gunzip -c "$BACKUP_FILE" | docker exec -i fidasim-postgres psql -U fidasim fidasim
    else
        docker exec -i fidasim-postgres psql -U fidasim fidasim < "$BACKUP_FILE"
    fi

    if [ $? -eq 0 ]; then
        print_info "Database restored successfully"
    else
        print_error "Database restore failed"
        exit 1
    fi
}

# Function to run database migrations
run_migrations() {
    print_info "Running database migrations..."

    # For now, just ensure the init script has run
    docker exec fidasim-postgres psql -U fidasim -d fidasim -c "SELECT COUNT(*) FROM users;" > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        print_info "Database schema is up to date"
    else
        print_warning "Running initial database setup..."
        docker exec -i fidasim-postgres psql -U fidasim fidasim < "$PROJECT_ROOT/docker/postgres/init.sql"
    fi
}

# Function to scale service
scale_service() {
    SERVICE="$1"
    REPLICAS="$2"

    if [ -z "$SERVICE" ] || [ -z "$REPLICAS" ]; then
        print_error "Usage: $0 scale <service> <replicas>"
        exit 1
    fi

    print_info "Scaling $SERVICE to $REPLICAS replicas..."

    cd "$PROJECT_ROOT"
    $COMPOSE_CMD up -d --scale "$SERVICE=$REPLICAS"

    print_info "Scaling complete"
}

# Parse command line arguments
case "$1" in
    start)
        ACTION="start"
        shift
        ;;
    stop)
        ACTION="stop"
        shift
        ;;
    restart)
        ACTION="restart"
        shift
        ;;
    status)
        ACTION="status"
        shift
        ;;
    logs)
        ACTION="logs"
        SERVICE="$2"
        shift 2
        ;;
    backup)
        ACTION="backup"
        shift
        ;;
    restore)
        ACTION="restore"
        BACKUP_FILE="$2"
        shift 2
        ;;
    migrate)
        ACTION="migrate"
        shift
        ;;
    scale)
        ACTION="scale"
        SERVICE="$2"
        REPLICAS="$3"
        shift 3
        ;;
    health)
        ACTION="health"
        shift
        ;;
    help|--help|-h)
        echo "Usage: $0 [COMMAND] [OPTIONS]"
        echo ""
        echo "Commands:"
        echo "  start              Start all services"
        echo "  stop               Stop all services"
        echo "  restart            Restart all services"
        echo "  status             Show service status"
        echo "  logs [service]     Show logs (all services or specific)"
        echo "  backup             Backup database"
        echo "  restore <file>     Restore database from backup"
        echo "  migrate            Run database migrations"
        echo "  scale <svc> <n>    Scale service to n replicas"
        echo "  health             Check service health"
        echo ""
        echo "Options:"
        echo "  --dev              Use development configuration"
        echo "  --prod             Use production configuration"
        echo ""
        echo "Examples:"
        echo "  $0 start --dev"
        echo "  $0 logs fidasim-core"
        echo "  $0 scale fidasim-core 3"
        echo "  $0 backup"
        exit 0
        ;;
    *)
        if [ -n "$1" ]; then
            print_error "Unknown command: $1"
            echo "Run '$0 help' for usage information"
            exit 1
        fi
        ACTION="start"
        ;;
esac

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --dev)
            DEPLOYMENT_MODE="development"
            shift
            ;;
        --prod)
            DEPLOYMENT_MODE="production"
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# Execute action
case "$ACTION" in
    start)
        check_environment
        start_services
        ;;
    stop)
        stop_services
        ;;
    restart)
        restart_services
        ;;
    status)
        show_status
        ;;
    logs)
        show_logs "$SERVICE"
        ;;
    backup)
        backup_database
        ;;
    restore)
        restore_database "$BACKUP_FILE"
        ;;
    migrate)
        run_migrations
        ;;
    scale)
        scale_service "$SERVICE" "$REPLICAS"
        ;;
    health)
        check_service_health
        ;;
esac