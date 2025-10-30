# FIDASIM Web Platform - Deployment Guide

Complete guide for deploying the FIDASIM containerized platform.

## Overview

The FIDASIM platform consists of 7 containerized services:

1. **fidasim-core** - Fortran simulation engine (Port 8001)
2. **fidasim-preprocessor** - Python preprocessing tools (Port 8002)
3. **fidasim-web-backend** - FastAPI REST API (Port 8000)
4. **fidasim-web-frontend** - React web interface (Port 80)
5. **postgres** - PostgreSQL database (Port 5433)
6. **redis** - Redis cache (Port 6379)

## Quick Start

### Prerequisites

- Docker 20.10+ with Docker Compose
- 8GB+ RAM
- 20GB+ free disk space
- Linux, macOS, or Windows with WSL2

### 1. Clone and Navigate

```bash
cd /path/to/FIDASIM
```

### 2. Configure Environment

```bash
# Copy environment template
cp .env.prod.example .env

# Generate a secure secret key
openssl rand -hex 32

# Edit .env and update:
# - SECRET_KEY (use the generated key above)
# - POSTGRES_PASSWORD (change from default)
# - VITE_API_URL and VITE_WS_URL (if deploying remotely)
```

**Important Environment Variables:**

```bash
# Security (REQUIRED)
SECRET_KEY=your-generated-secret-key-here
POSTGRES_PASSWORD=change-this-secure-password

# URLs (update for your domain)
VITE_API_URL=http://your-domain.com:8000  # or https://your-domain.com/api
VITE_WS_URL=ws://your-domain.com:8000     # or wss://your-domain.com/ws
```

### 3. Build and Deploy

```bash
# One-command deployment (builds and starts everything)
./scripts/deploy-web.sh
```

Or step by step:

```bash
# Build images
./scripts/build-web.sh

# Start all services
docker-compose -f docker-compose.prod.yaml up -d

# Test deployment
./scripts/test-web.sh
```

### 4. Access the Platform

Once deployed, access:

- **Web Interface**: http://localhost (or your domain)
- **Backend API**: http://localhost:8000
- **API Documentation**: http://localhost:8000/api/docs
- **Core API**: http://localhost:8001/docs
- **Preprocessor API**: http://localhost:8002/docs

### 5. Create Initial Admin User

The system requires an initial user. Two options:

**Option A: Using Docker exec**

```bash
docker exec -it fidasim-web-backend python -c "
from app.core.security import get_password_hash
from app.models.user import User
from app.db.database import SessionLocal

db = SessionLocal()
admin = User(
    username='admin',
    password_hash=get_password_hash('changeme'),
    role='admin',
    email='admin@example.com',
    is_active=True
)
db.add(admin)
db.commit()
print('Admin user created!')
"
```

**Option B: Using SQL**

```bash
docker exec -it fidasim-postgres psql -U fidasim -d fidasim -c "
INSERT INTO users (id, username, password_hash, role, email, is_active, max_cores, max_jobs)
VALUES (
    gen_random_uuid(),
    'admin',
    '\$2b\$12\$...',  -- Use bcrypt hash of your password
    'admin',
    'admin@example.com',
    true,
    32,
    10
);
"
```

Then login at http://localhost with:
- Username: `admin`
- Password: `changeme` (or your chosen password)

## Architecture

```
┌─────────────────────────────────────────────────┐
│  Browser (User)                                 │
└────────────────┬────────────────────────────────┘
                 │ HTTP/WebSocket
                 v
┌─────────────────────────────────────────────────┐
│  Nginx (fidasim-web-frontend)                   │
│  Port 80 - Serves React app, proxies API       │
└────────┬────────────────────────────────────────┘
         │
         ├── /api/* ──────────────┐
         │                        v
         │               ┌────────────────────────┐
         │               │  FastAPI Backend       │
         │               │  Port 8000             │
         │               │  - REST API            │
         │               │  - WebSocket           │
         │               │  - Job orchestration   │
         │               └───┬───────────┬────────┘
         │                   │           │
         v                   v           v
    ┌─────────┐      ┌──────────┐  ┌─────────┐
    │ Core    │      │PostgreSQL│  │  Redis  │
    │ Port    │      │Port 5433 │  │Port 6379│
    │ 8001    │      └──────────┘  └─────────┘
    └─────────┘
         │
    ┌─────────┐
    │Preproc  │
    │Port 8002│
    └─────────┘
```

## Service Details

### Backend (fidasim-web-backend)

**Image**: `fidasim-web-backend:latest`
**Base**: Python 3.10-slim
**Workers**: 4 (uvicorn)
**Health Check**: `GET /api/health`

**Environment Variables**:
- `SECRET_KEY` - JWT signing key (REQUIRED)
- `DATABASE_URL` - PostgreSQL connection string
- `REDIS_URL` - Redis connection string
- `CORE_API_URL` - URL to fidasim-core
- `PREPROCESSOR_API_URL` - URL to preprocessor
- `MAX_CORES`, `MAX_UPLOAD_SIZE` - Resource limits

**Volumes**:
- `/data` - Shared simulation data
- `/data/uploads` - User-uploaded files

### Frontend (fidasim-web-frontend)

**Image**: `fidasim-web-frontend:latest`
**Base**: Nginx Alpine
**Multi-stage**: Node 18 (build) → Nginx (serve)

**Build Args** (compile-time):
- `VITE_API_URL`
- `VITE_WS_URL`
- `VITE_APP_NAME`
- `VITE_APP_VERSION`

**Environment Variables** (runtime):
- Same as build args, injected at container start

**Features**:
- Static file caching
- API proxy to backend
- WebSocket upgrade support
- Gzip compression
- Security headers

## Management Commands

### View Logs

```bash
# All services
docker-compose -f docker-compose.prod.yaml logs -f

# Specific service
docker logs -f fidasim-web-backend
docker logs -f fidasim-web-frontend

# Last 100 lines
docker logs --tail=100 fidasim-web-backend
```

### Service Control

```bash
# Stop all services
docker-compose -f docker-compose.prod.yaml down

# Stop web services only
docker-compose -f docker-compose.prod.yaml stop fidasim-web-backend fidasim-web-frontend

# Restart a service
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend

# Rebuild and restart
docker-compose -f docker-compose.prod.yaml up -d --build fidasim-web-backend
```

### Database Operations

```bash
# Access PostgreSQL
docker exec -it fidasim-postgres psql -U fidasim -d fidasim

# Backup database
docker exec fidasim-postgres pg_dump -U fidasim fidasim > backup.sql

# Restore database
cat backup.sql | docker exec -i fidasim-postgres psql -U fidasim -d fidasim
```

### Redis Operations

```bash
# Access Redis CLI
docker exec -it fidasim-redis redis-cli

# View keys
docker exec fidasim-redis redis-cli KEYS '*'

# Flush cache
docker exec fidasim-redis redis-cli FLUSHALL
```

## Troubleshooting

### Backend Won't Start

**Check logs:**
```bash
docker logs fidasim-web-backend
```

**Common issues:**
1. Database not ready - Wait 30s and restart
2. SECRET_KEY not set - Check .env file
3. Port 8000 in use - Stop conflicting service

**Solution:**
```bash
# Wait for database
docker-compose -f docker-compose.prod.yaml up -d postgres
sleep 30

# Restart backend
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend
```

### Frontend Won't Load

**Check logs:**
```bash
docker logs fidasim-web-frontend
```

**Common issues:**
1. Backend not accessible - Check VITE_API_URL
2. Build failed - Check build logs
3. Nginx config error - Check nginx.conf syntax

**Solution:**
```bash
# Rebuild frontend
docker-compose -f docker-compose.prod.yaml build --no-cache fidasim-web-frontend
docker-compose -f docker-compose.prod.yaml up -d fidasim-web-frontend
```

### WebSocket Connection Fails

**Check:**
1. VITE_WS_URL matches backend URL (use `ws://` not `http://`)
2. Nginx WebSocket proxy configured correctly
3. Firewall allows WebSocket connections

**Test WebSocket:**
```bash
# From inside backend container
docker exec -it fidasim-web-backend python -c "
from fastapi.testclient import TestClient
from app.main import app
client = TestClient(app)
# Test WebSocket endpoint
"
```

### Database Connection Issues

**Check connection:**
```bash
docker exec fidasim-web-backend python -c "
from app.db.database import engine
try:
    engine.connect()
    print('Database connection OK')
except Exception as e:
    print(f'Database connection failed: {e}')
"
```

**Reset database:**
```bash
# WARNING: This deletes all data!
docker-compose -f docker-compose.prod.yaml down -v
docker-compose -f docker-compose.prod.yaml up -d postgres
# Wait for init
sleep 30
docker-compose -f docker-compose.prod.yaml up -d
```

## Production Deployment

### Security Checklist

- [ ] Change `SECRET_KEY` to a securely generated value
- [ ] Change `POSTGRES_PASSWORD` from default
- [ ] Update CORS origins in backend config
- [ ] Use HTTPS/WSS in production (not HTTP/WS)
- [ ] Set up SSL certificates
- [ ] Configure firewall rules
- [ ] Disable debug mode
- [ ] Set resource limits in docker-compose
- [ ] Regular database backups
- [ ] Monitor logs and alerts

### Using HTTPS/WSS

1. **Get SSL certificate** (Let's Encrypt recommended)

2. **Update Nginx configuration**:
```nginx
server {
    listen 443 ssl http2;
    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;
    # ... rest of config
}
```

3. **Update environment variables**:
```bash
VITE_API_URL=https://your-domain.com/api
VITE_WS_URL=wss://your-domain.com/api/ws
```

4. **Rebuild frontend** with new URLs

### Resource Limits

Add to docker-compose.prod.yaml:

```yaml
services:
  fidasim-web-backend:
    deploy:
      resources:
        limits:
          cpus: '4'
          memory: 4G
        reservations:
          cpus: '2'
          memory: 2G
```

### Monitoring

**Health checks:**
```bash
# Automated health check script
watch -n 30 './scripts/test-web.sh'
```

**Prometheus metrics** (optional):
- Add prometheus exporters
- Configure Grafana dashboards
- Set up alerts

## Backup and Recovery

### Backup Script

```bash
#!/bin/bash
DATE=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="/backups/$DATE"

mkdir -p "$BACKUP_DIR"

# Backup database
docker exec fidasim-postgres pg_dump -U fidasim fidasim > "$BACKUP_DIR/database.sql"

# Backup volumes
docker run --rm -v fidasim-data:/data -v "$BACKUP_DIR:/backup" alpine tar czf /backup/data.tar.gz /data
docker run --rm -v web-uploads:/uploads -v "$BACKUP_DIR:/backup" alpine tar czf /backup/uploads.tar.gz /uploads

echo "Backup completed: $BACKUP_DIR"
```

### Recovery

```bash
# Restore database
cat backup/database.sql | docker exec -i fidasim-postgres psql -U fidasim -d fidasim

# Restore volumes
docker run --rm -v fidasim-data:/data -v "$BACKUP_DIR:/backup" alpine tar xzf /backup/data.tar.gz -C /
```

## Updating

### Update Backend Code

```bash
cd docker/web/backend
# Make code changes
docker-compose -f ../../docker-compose.prod.yaml build fidasim-web-backend
docker-compose -f ../../docker-compose.prod.yaml up -d fidasim-web-backend
```

### Update Frontend Code

```bash
cd docker/web/frontend
# Make code changes
docker-compose -f ../../docker-compose.prod.yaml build fidasim-web-frontend
docker-compose -f ../../docker-compose.prod.yaml up -d fidasim-web-frontend
```

## Support

For issues:
1. Check logs: `docker-compose logs`
2. Run tests: `./scripts/test-web.sh`
3. Check this guide's troubleshooting section
4. Review container health: `docker ps`

## Performance Tuning

### Backend
- Increase uvicorn workers (default: 4)
- Adjust `MAX_CONCURRENT_JOBS`
- Configure PostgreSQL connection pool
- Optimize Redis memory settings

### Frontend
- Enable Nginx caching
- Configure CDN for static assets
- Optimize image sizes
- Use HTTP/2

### Database
- Add indexes for frequent queries
- Configure PostgreSQL shared_buffers
- Set up read replicas for scaling
- Regular VACUUM operations

---

**Last Updated**: Phase 2 Completion
**Version**: 1.0.0
