# FIDASIM Containerized Application - Phase 1 Documentation

## Overview

This document describes the containerized FIDASIM application architecture and provides instructions for building, deploying, and using the system. Phase 1 implements the core containers for running FIDASIM simulations and preprocessing data.

## Architecture

The system consists of the following containers:

### Core Services (Phase 1 - Implemented)
- **fidasim-core**: Fortran simulation engine with REST API
- **fidasim-preprocessor**: Python preprocessing service for data preparation
- **postgres**: PostgreSQL database for job tracking and metadata
- **redis**: Redis for caching and job queue management

### Web Interface (Phase 2+ - To Be Implemented)
- **fidasim-web**: Web interface for job submission and monitoring

## Prerequisites

- Docker 20.10+ and Docker Compose 2.0+
- At least 16GB RAM and 50GB disk space
- Linux or macOS (WSL2 for Windows)
- Git for cloning the repository

## Quick Start

### 1. Clone and Setup

```bash
# Clone the repository
git clone <repository-url>
cd FIDASIM

# Copy environment template
cp .env.template .env

# Edit .env and update passwords
nano .env
```

### 2. Build Docker Images

```bash
# Build all images (first build takes ~45-60 minutes)
./scripts/build.sh

# Or build with specific options
./scripts/build.sh --dev          # Development mode
./scripts/build.sh --no-cache     # Force rebuild
./scripts/build.sh --test         # Run tests after build
```

### 3. Deploy Services

```bash
# Start all services
./scripts/deploy.sh start

# Or start in development mode
./scripts/deploy.sh start --dev

# Check service health
./scripts/deploy.sh health
```

### 4. Verify Installation

```bash
# Check service status
./scripts/deploy.sh status

# View logs
./scripts/deploy.sh logs

# Test API endpoints
curl http://localhost:8001/health  # Core API
curl http://localhost:8002/health  # Preprocessor API
```

## API Documentation

### FIDASIM Core API (Port 8001)

The core API provides endpoints for running FIDASIM simulations:

#### Endpoints

**Health Check**
```bash
GET /health
```

**Run Simulation**
```bash
POST /run
Content-Type: application/json

{
  "runid": "test_run_001",
  "input_file": "path/to/namelist.dat",
  "cores": 8,
  "equilibrium_file": "path/to/equilibrium.h5",
  "geometry_file": "path/to/geometry.h5",
  "distribution_file": "path/to/distribution.h5"
}
```

**Check Job Status**
```bash
GET /status/{job_id}
```

**Get Job Logs**
```bash
GET /logs/{job_id}?tail=100&stream=false
```

**List Jobs**
```bash
GET /jobs?status=running&limit=100
```

**Download Results**
```bash
GET /results/{job_id}/{filename}
```

**Upload File**
```bash
POST /upload
Content-Type: multipart/form-data

file: <file>
directory: <optional-subdirectory>
```

### Preprocessor API (Port 8002)

The preprocessor API provides data preparation services:

#### Endpoints

**Run PREFIDA**
```bash
POST /prefida
Content-Type: multipart/form-data

config: {
  "shot": 100000,
  "time": 1.0,
  "device": "DIII-D",
  "runid": "prefida_001",
  "nr": 100,
  "nz": 100,
  "calc_bes": true,
  "calc_fida": true
}
equilibrium_file: <GEQDSK file>
profiles_file: <optional>
nbi_file: <optional>
```

**Convert GEQDSK**
```bash
POST /convert/geqdsk
Content-Type: multipart/form-data

file: <GEQDSK file>
output_format: "hdf5"
```

**Generate Grid**
```bash
POST /grid/generate
Content-Type: application/json

{
  "grid_type": "rz_uniform",
  "nr": 100,
  "nz": 100,
  "nphi": 1,
  "rmin": 0.5,
  "rmax": 2.5,
  "zmin": -1.5,
  "zmax": 1.5
}
```

**Get Templates**
```bash
GET /templates/namelist
GET /templates/prefida_config
GET /templates/grid_config
```

## Usage Examples

### Example 1: Running a Basic FIDASIM Simulation

```python
import requests
import json

# Upload input files
with open('namelist.dat', 'rb') as f:
    response = requests.post(
        'http://localhost:8001/upload',
        files={'file': f},
        data={'directory': 'run_001'}
    )

# Start simulation
config = {
    "runid": "test_001",
    "input_file": "run_001/namelist.dat",
    "cores": 8
}

response = requests.post(
    'http://localhost:8001/run',
    json=config
)

job = response.json()
job_id = job['job_id']

# Check status
response = requests.get(f'http://localhost:8001/status/{job_id}')
print(response.json())

# Stream logs
response = requests.get(
    f'http://localhost:8001/logs/{job_id}',
    params={'stream': True}
)
for line in response.iter_lines():
    print(line.decode())
```

### Example 2: Preprocessing with PREFIDA

```python
import requests

# Prepare files
files = {
    'equilibrium_file': open('g000001.01000', 'rb'),
    'profiles_file': open('profiles.dat', 'rb')
}

# Configure PREFIDA
config = {
    "shot": 100000,
    "time": 1.0,
    "device": "DIII-D",
    "runid": "prefida_test",
    "nr": 100,
    "nz": 100,
    "calc_bes": True,
    "calc_fida": True
}

# Run preprocessing
response = requests.post(
    'http://localhost:8002/prefida',
    data={'config': json.dumps(config)},
    files=files
)

result = response.json()
if result['success']:
    print(f"Output file: {result['output_file']}")
    print(f"Metadata: {result['metadata']}")
```

### Example 3: Using Docker Compose Directly

```bash
# Start specific service
docker-compose up -d fidasim-core

# Scale service
docker-compose up -d --scale fidasim-core=3

# Execute command in container
docker exec fidasim-core /app/fidasim /data/test.dat

# View real-time logs
docker-compose logs -f fidasim-core

# Stop services
docker-compose down
```

## Directory Structure

```
FIDASIM/
├── docker/
│   ├── core/
│   │   ├── Dockerfile          # Core service Docker image
│   │   ├── api.py             # FastAPI application
│   │   └── entrypoint.sh      # Container entrypoint
│   ├── preprocessor/
│   │   ├── Dockerfile          # Preprocessor Docker image
│   │   ├── api.py             # FastAPI application
│   │   └── requirements.txt   # Python dependencies
│   └── postgres/
│       └── init.sql           # Database schema
├── data/                      # Shared data volume
├── scripts/
│   ├── build.sh              # Build script
│   └── deploy.sh             # Deployment script
├── docker-compose.yaml        # Main compose file
├── docker-compose.dev.yaml    # Development overrides
└── .env.template             # Environment template
```

## Database Schema

The PostgreSQL database tracks jobs, users, and performance metrics:

### Main Tables
- **users**: User accounts and permissions
- **jobs**: Job tracking and metadata
- **parameter_scans**: Parameter scan configurations
- **job_logs**: Detailed job execution logs
- **performance_metrics**: Performance tracking
- **files**: File metadata

### Default Credentials
- Database: `fidasim` / `fidasim`
- Admin user: `admin` / `admin123` (CHANGE THIS!)

## Development

### Building for Development

```bash
# Build with development settings
./scripts/build.sh --dev

# Start with hot-reload enabled
./scripts/deploy.sh start --dev

# Access development tools
http://localhost:8080  # Adminer (database UI)
http://localhost:8081  # Redis Commander
```

### Running Tests

```bash
# Test containers after build
./scripts/build.sh --test

# Run specific container tests
docker run --rm fidasim-core:latest /app/entrypoint.sh test

# Test API endpoints
curl -X GET http://localhost:8001/health
curl -X GET http://localhost:8002/health
```

### Debugging

```bash
# Enter container shell
docker exec -it fidasim-core bash

# View detailed logs
docker logs -f fidasim-core --tail 100

# Check resource usage
docker stats fidasim-core

# Inspect container
docker inspect fidasim-core
```

## Troubleshooting

### Common Issues

**1. Build Fails - Atomic Tables**
```bash
# Pre-build atomic tables locally
cd /path/to/FIDASIM
make tables  # Takes 30+ minutes
```

**2. Port Already in Use**
```bash
# Check what's using the port
lsof -i :8001

# Change port in docker-compose.yaml
ports:
  - "8091:8001"  # Map to different host port
```

**3. Permission Denied**
```bash
# Fix data directory permissions
chmod 755 data/
chown -R $(id -u):$(id -g) data/
```

**4. Out of Memory**
```bash
# Increase Docker memory limit
# Docker Desktop: Settings > Resources > Memory
# Or in docker-compose.yaml:
deploy:
  resources:
    limits:
      memory: 16G
```

**5. Database Connection Failed**
```bash
# Reset database
docker-compose down -v
docker-compose up -d postgres
./scripts/deploy.sh migrate
```

## Performance Tuning

### Docker Settings
```yaml
# In docker-compose.yaml
deploy:
  resources:
    limits:
      cpus: '32'
      memory: 8G
    reservations:
      cpus: '4'
      memory: 2G
```

### OpenMP Threads
```bash
# Set in .env
MAX_CORES=32
DEFAULT_CORES_PER_JOB=8

# Or per job
curl -X POST http://localhost:8001/run \
  -H "Content-Type: application/json" \
  -d '{"cores": 16, ...}'
```

### Database Optimization
```sql
-- Adjust in postgres/init.sql
-- Connection pool size
DB_POOL_SIZE=20

-- Work memory
ALTER SYSTEM SET work_mem = '256MB';
ALTER SYSTEM SET shared_buffers = '2GB';
```

## Security Considerations

### Production Deployment

1. **Change Default Passwords**
```bash
# Generate secure passwords
openssl rand -base64 32  # For database
openssl rand -hex 32      # For secret keys
```

2. **Enable HTTPS**
```nginx
# Add nginx reverse proxy with SSL
server {
    listen 443 ssl;
    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;

    location / {
        proxy_pass http://localhost:8000;
    }
}
```

3. **Network Isolation**
```yaml
# Use internal networks
networks:
  internal:
    internal: true
  external:
    internal: false
```

4. **Resource Limits**
```yaml
# Enforce limits
deploy:
  resources:
    limits:
      cpus: '4'
      memory: 4G
```

## Monitoring

### Health Checks
```bash
# Automated health monitoring
watch -n 5 './scripts/deploy.sh health'
```

### Metrics Collection
```bash
# Prometheus endpoint (when implemented)
curl http://localhost:9090/metrics
```

### Log Aggregation
```bash
# Collect all logs
docker-compose logs --tail=100 > fidasim.log
```

## Backup and Recovery

### Backup Database
```bash
./scripts/deploy.sh backup
```

### Restore Database
```bash
./scripts/deploy.sh restore backups/fidasim_db_20240101_120000.sql.gz
```

### Backup Data Files
```bash
tar -czf fidasim_data_$(date +%Y%m%d).tar.gz data/
```

## Next Steps (Phase 2+)

The following features will be implemented in subsequent phases:

1. **Web Interface**
   - React frontend for job submission
   - Real-time job monitoring dashboard
   - Parameter scan builder
   - Result visualization

2. **Advanced Features**
   - Smart core allocation algorithm
   - WebSocket support for live updates
   - OAuth2/LDAP authentication
   - S3/object storage integration

3. **Kubernetes Migration**
   - Helm charts
   - Horizontal pod autoscaling
   - Ingress configuration
   - Persistent volume claims

## Support

For issues or questions:
- Check the troubleshooting section above
- Review logs: `./scripts/deploy.sh logs`
- File an issue on GitHub
- Contact the development team

## License

See LICENSE.md in the repository root.