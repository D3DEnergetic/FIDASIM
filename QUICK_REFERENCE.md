# FIDASIM Platform - Quick Reference Card

## 🚀 Quick Start (3 Steps)

```bash
# 1. Configure (first time only)
cp .env.prod.example .env
# Edit .env: Set SECRET_KEY and POSTGRES_PASSWORD

# 2. Deploy everything
./scripts/deploy-web.sh

# 3. Access
open http://localhost
```

## 📦 Management Commands

### Build & Deploy
```bash
# Build images
./scripts/build-web.sh

# Deploy all services
./scripts/deploy-web.sh

# Test deployment
./scripts/test-web.sh
```

### Service Control
```bash
# Start all
docker-compose -f docker-compose.prod.yaml up -d

# Stop all
docker-compose -f docker-compose.prod.yaml down

# Restart service
docker-compose -f docker-compose.prod.yaml restart <service>

# View status
docker-compose -f docker-compose.prod.yaml ps
```

### Logs
```bash
# All services
docker-compose -f docker-compose.prod.yaml logs -f

# Specific service
docker logs -f fidasim-web-backend
docker logs -f fidasim-web-frontend

# Last 100 lines
docker logs --tail=100 fidasim-web-backend
```

## 🔧 Common Tasks

### Create Admin User
```bash
docker exec -it fidasim-web-backend python -c "
from app.core.security import get_password_hash
from app.models.user import User
from app.db.database import SessionLocal
db = SessionLocal()
admin = User(
    username='admin',
    password_hash=get_password_hash('admin123'),
    role='admin',
    email='admin@example.com',
    is_active=True
)
db.add(admin)
db.commit()
print('✓ Admin created: admin/admin123')
"
```

### Database Operations
```bash
# Access PostgreSQL
docker exec -it fidasim-postgres psql -U fidasim -d fidasim

# Backup
docker exec fidasim-postgres pg_dump -U fidasim fidasim > backup.sql

# Restore
cat backup.sql | docker exec -i fidasim-postgres psql -U fidasim -d fidasim
```

### Redis Operations
```bash
# Access Redis
docker exec -it fidasim-redis redis-cli

# Check connection
docker exec fidasim-redis redis-cli PING

# Flush cache
docker exec fidasim-redis redis-cli FLUSHALL
```

## 🌐 Access Points

| Service | URL | Description |
|---------|-----|-------------|
| **Web UI** | http://localhost | Main interface |
| **Backend API** | http://localhost:8000 | REST API |
| **API Docs** | http://localhost:8000/api/docs | Swagger UI |
| **Core API** | http://localhost:8001/docs | Simulation API |
| **Preprocessor** | http://localhost:8002/docs | Preprocessing API |

## 🐛 Troubleshooting

### Backend Won't Start
```bash
# Check logs
docker logs fidasim-web-backend

# Common fix: Restart with fresh DB connection
docker-compose -f docker-compose.prod.yaml restart postgres
sleep 30
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend
```

### Frontend Won't Load
```bash
# Check logs
docker logs fidasim-web-frontend

# Rebuild
docker-compose -f docker-compose.prod.yaml build --no-cache fidasim-web-frontend
docker-compose -f docker-compose.prod.yaml up -d fidasim-web-frontend
```

### WebSocket Not Connecting
```bash
# Check backend is running
curl http://localhost:8000/api/health

# Check frontend can reach backend
docker exec fidasim-web-frontend wget -q -O- http://fidasim-web-backend:8000/api/health
```

### Database Connection Failed
```bash
# Check database is ready
docker exec fidasim-postgres pg_isready -U fidasim

# If not ready, wait and restart
sleep 30
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend
```

## 🔐 Environment Variables

### Required
```bash
SECRET_KEY=<generate with: openssl rand -hex 32>
POSTGRES_PASSWORD=<secure password>
```

### Optional (with defaults)
```bash
VITE_API_URL=http://localhost:8000
VITE_WS_URL=ws://localhost:8000
MAX_CORES=32
MAX_UPLOAD_SIZE=1073741824
```

## 📊 Service Status

### Check Health
```bash
# All services
./scripts/test-web.sh

# Individual health checks
curl http://localhost:8000/api/health  # Backend
curl http://localhost:8001/health       # Core
curl http://localhost:8002/health       # Preprocessor
curl http://localhost/                  # Frontend
```

### View Resources
```bash
# Container stats
docker stats

# Disk usage
docker system df

# Network info
docker network inspect fidasim-network
```

## 🔄 Updates

### Update Backend
```bash
cd docker/web/backend
# Make changes to code
docker-compose -f ../../docker-compose.prod.yaml build fidasim-web-backend
docker-compose -f ../../docker-compose.prod.yaml up -d fidasim-web-backend
```

### Update Frontend
```bash
cd docker/web/frontend
# Make changes to code
docker-compose -f ../../docker-compose.prod.yaml build fidasim-web-frontend
docker-compose -f ../../docker-compose.prod.yaml up -d fidasim-web-frontend
```

### Update Configuration
```bash
# Edit .env file
nano .env

# Restart services to apply
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend fidasim-web-frontend
```

## 🧹 Cleanup

### Remove Stopped Containers
```bash
docker-compose -f docker-compose.prod.yaml down
```

### Remove Volumes (⚠️  Deletes data!)
```bash
docker-compose -f docker-compose.prod.yaml down -v
```

### Clean Docker System
```bash
# Remove unused images, containers, networks
docker system prune -a

# Remove unused volumes
docker volume prune
```

## 📝 Useful Scripts

| Script | Purpose |
|--------|---------|
| `./scripts/build-web.sh` | Build web service images |
| `./scripts/deploy-web.sh` | Deploy complete stack |
| `./scripts/test-web.sh` | Test all services |

## 📚 Documentation

| File | Description |
|------|-------------|
| `PROJECT_COMPLETE.md` | Full project summary |
| `DEPLOYMENT_GUIDE.md` | Detailed deployment guide |
| `DOCKERIZATION_COMPLETE.md` | Dockerization details |
| `FRONTEND_COMPLETE.md` | Frontend documentation |
| `BACKEND_COMPLETE.md` | Backend documentation |

## 🆘 Getting Help

1. Check logs: `docker logs <container>`
2. Run tests: `./scripts/test-web.sh`
3. Review docs: `DEPLOYMENT_GUIDE.md`
4. Check status: `docker ps`
5. Inspect network: `docker network inspect fidasim-network`

## 🎯 Common Workflows

### Deploy Fresh Installation
```bash
cp .env.prod.example .env
# Edit .env
./scripts/deploy-web.sh
# Create admin user (see above)
# Login at http://localhost
```

### Daily Operations
```bash
# Check status
docker-compose -f docker-compose.prod.yaml ps

# View logs
docker-compose -f docker-compose.prod.yaml logs -f --tail=50

# Monitor resources
docker stats
```

### Restart Everything
```bash
docker-compose -f docker-compose.prod.yaml restart
```

### Backup Before Changes
```bash
# Backup database
docker exec fidasim-postgres pg_dump -U fidasim fidasim > backup_$(date +%Y%m%d).sql

# Backup volumes
docker run --rm -v fidasim-data:/data -v $(pwd):/backup alpine tar czf /backup/data_$(date +%Y%m%d).tar.gz /data
```

## ⚡ Performance Tips

- **Backend**: Increase workers in Dockerfile (default: 4)
- **Frontend**: Enable CDN for static assets
- **Database**: Tune PostgreSQL settings for your workload
- **Redis**: Increase memory if needed
- **Nginx**: Enable caching, compression (already configured)

## 🔒 Security Checklist

- [ ] Changed `SECRET_KEY` from default
- [ ] Changed `POSTGRES_PASSWORD` from default
- [ ] Using HTTPS/WSS in production
- [ ] Firewall configured
- [ ] Regular backups enabled
- [ ] Monitoring alerts set up
- [ ] Log rotation configured
- [ ] Updated CORS origins

---

**Quick Deploy**: `./scripts/deploy-web.sh`
**Access**: http://localhost
**Docs**: http://localhost:8000/api/docs
**Status**: `docker-compose -f docker-compose.prod.yaml ps`
**Logs**: `docker-compose -f docker-compose.prod.yaml logs -f`
