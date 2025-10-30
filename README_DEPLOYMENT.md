# FIDASIM Web Platform - Deployment README

Complete web interface for FIDASIM (Fast-Ion Diagnostic Simulation) with job submission, monitoring, and interactive data visualization.

## Features

- **Web-Based Job Submission** - Submit FIDASIM simulations through a user-friendly interface
- **Real-Time Monitoring** - Watch jobs progress with live status updates via WebSocket
- **Interactive Visualization** - Explore multi-dimensional HDF5/NetCDF output with Plotly
  - 1D line plots and 2D heatmaps
  - Dimension slicing and integration
  - Dynamic coordinate selection
  - Export capabilities
- **User Management** - Role-based access control (admin/user)
- **File Management** - Upload input files, download results
- **Job History** - Track all simulations with searchable history

## Quick Start

### Prerequisites

- Docker 20.10+ with Docker Compose
- 8GB+ RAM
- 20GB+ free disk space
- Linux, macOS, or Windows with WSL2

### 1. Clone Repository

```bash
git clone https://github.com/yourusername/FIDASIM.git
cd FIDASIM
```

### 2. Configure Environment

```bash
# Copy environment template
cp .env.prod.example .env

# Generate secure keys
openssl rand -hex 32  # Use this for SECRET_KEY

# Edit .env and set:
# - SECRET_KEY=<generated key>
# - POSTGRES_PASSWORD=<secure password>
# - VITE_API_URL=http://localhost:8000 (or your domain)
# - VITE_WS_URL=ws://localhost:8000 (or your domain)
```

### 3. Deploy

```bash
# One-command deployment
./scripts/deploy-web.sh

# Or step by step:
./scripts/build-web.sh
docker-compose -f docker-compose.prod.yaml up -d
./scripts/test-web.sh
```

### 4. Create Admin User

```bash
docker exec fidasim-web-backend python -c "
from app.core.security import get_password_hash
from app.models.user import User
from app.db.database import SessionLocal

db = SessionLocal()
admin = User(
    username='admin',
    password_hash=get_password_hash('admin123'),
    role='admin',
    email='admin@example.com',
    is_active=True,
    max_cores=32,
    max_jobs=10
)
db.add(admin)
db.commit()
print('✓ Admin user created: admin/admin123')
"
```

### 5. Access Platform

- **Web Interface**: http://localhost:3000
- **Backend API**: http://localhost:8000
- **API Docs**: http://localhost:8000/api/docs

Login with:
- Username: `admin`
- Password: `admin123` (change this!)

## Architecture

```
┌─────────────────────────────────────────────────┐
│  Browser                                        │
└────────────────┬────────────────────────────────┘
                 │
                 v
┌─────────────────────────────────────────────────┐
│  Nginx (Frontend) - Port 3000                   │
│  Serves React app + API proxy                   │
└────────┬────────────────────────────────────────┘
         │
         v
┌─────────────────────────────────────────────────┐
│  FastAPI Backend - Port 8000                    │
│  REST API + WebSocket + Job orchestration       │
└───┬─────────┬─────────┬─────────────────────────┘
    │         │         │
    v         v         v
┌─────┐  ┌─────┐  ┌────────┐
│Core │  │Pre- │  │Postgres│
│8001 │  │proc │  │  5433  │
│     │  │8002 │  │        │
└─────┘  └─────┘  └────────┘
```

## Services

1. **fidasim-web-frontend** (Port 3000) - React UI with Nginx
2. **fidasim-web-backend** (Port 8000) - FastAPI REST API
3. **fidasim-core** (Port 8001) - FIDASIM simulation engine
4. **fidasim-preprocessor** (Port 8002) - Input preprocessing
5. **postgres** (Port 5433) - PostgreSQL database
6. **redis** (Port 6379) - Redis cache

## Common Commands

```bash
# View logs
docker-compose -f docker-compose.prod.yaml logs -f

# Restart services
docker-compose -f docker-compose.prod.yaml restart

# Stop all
docker-compose -f docker-compose.prod.yaml down

# Rebuild after code changes
docker-compose -f docker-compose.prod.yaml build fidasim-web-backend
docker-compose -f docker-compose.prod.yaml up -d

# Check status
docker-compose -f docker-compose.prod.yaml ps
```

## Directory Structure

```
FIDASIM/
├── docker/
│   ├── web/
│   │   ├── backend/          # FastAPI application
│   │   │   ├── app/
│   │   │   │   ├── api/      # REST endpoints
│   │   │   │   ├── services/ # Business logic (incl. visualization)
│   │   │   │   ├── models/   # Database models
│   │   │   │   └── core/     # Config, security
│   │   │   ├── Dockerfile
│   │   │   └── requirements.txt
│   │   └── frontend/         # React application
│   │       ├── src/
│   │       │   ├── components/
│   │       │   ├── pages/
│   │       │   └── services/
│   │       ├── Dockerfile
│   │       └── package.json
│   ├── core/                 # FIDASIM simulation service
│   │   ├── Dockerfile        # (Python deps installed inline)
│   │   ├── api.py
│   │   └── entrypoint.sh
│   ├── preprocessor/         # Preprocessing service
│   │   ├── Dockerfile
│   │   ├── api.py
│   │   └── requirements.txt
│   └── postgres/             # Database
├── scripts/
│   ├── deploy-web.sh
│   ├── build-web.sh
│   └── test-web.sh
├── docker-compose.prod.yaml
├── .env.prod.example
└── DEPLOYMENT_GUIDE.md
```

## Using the Platform

### Submit a Job

1. Login to http://localhost:3000
2. Click "Submit Job"
3. Fill in job details:
   - Run ID (unique identifier)
   - Number of cores
   - Input file (upload or select existing)
4. Click "Submit"

### Monitor Jobs

1. Navigate to "My Jobs"
2. View real-time status updates
3. Click job to see details and logs
4. When complete, click "Visualize Results"

### Visualize Data

1. Select output file (HDF5/NetCDF)
2. Choose variable to plot
3. Select dimensions for axes
4. For other dimensions:
   - **Slice**: Use slider to select specific index
   - **Integrate**: Check box to sum over entire dimension
5. Click "Visualize Data"
6. Interact with plot (zoom, pan, export)

## Production Deployment

### Security Checklist

- [ ] Change `SECRET_KEY` to randomly generated value
- [ ] Change `POSTGRES_PASSWORD` from default
- [ ] Change admin password from `admin123`
- [ ] Use HTTPS/WSS (not HTTP/WS)
- [ ] Configure firewall rules
- [ ] Set up SSL certificates (Let's Encrypt)
- [ ] Update CORS origins in backend
- [ ] Enable regular database backups
- [ ] Set up monitoring and alerts

### HTTPS Setup

1. Obtain SSL certificate (Let's Encrypt recommended)
2. Update Nginx configuration with SSL
3. Update environment variables:
   ```bash
   VITE_API_URL=https://yourdomain.com/api
   VITE_WS_URL=wss://yourdomain.com/api/ws
   ```
4. Rebuild frontend with new URLs

### Resource Limits

Add to `docker-compose.prod.yaml`:

```yaml
services:
  fidasim-web-backend:
    deploy:
      resources:
        limits:
          cpus: '4'
          memory: 4G
```

## Troubleshooting

### Backend Won't Start

```bash
# Check logs
docker logs fidasim-web-backend

# Common fix: restart database
docker-compose -f docker-compose.prod.yaml restart postgres
sleep 30
docker-compose -f docker-compose.prod.yaml restart fidasim-web-backend
```

### Frontend Won't Load

```bash
# Check logs
docker logs fidasim-web-frontend

# Rebuild if needed
docker-compose -f docker-compose.prod.yaml build --no-cache fidasim-web-frontend
docker-compose -f docker-compose.prod.yaml up -d fidasim-web-frontend
```

### Database Issues

```bash
# Check connection
docker exec fidasim-postgres pg_isready -U fidasim

# Reset database (WARNING: deletes all data!)
docker-compose -f docker-compose.prod.yaml down -v
docker-compose -f docker-compose.prod.yaml up -d
```

## Backup and Recovery

### Backup Database

```bash
docker exec fidasim-postgres pg_dump -U fidasim fidasim > backup_$(date +%Y%m%d).sql
```

### Backup Job Data

```bash
docker run --rm -v fidasim-data:/data -v $(pwd):/backup alpine \
  tar czf /backup/data_$(date +%Y%m%d).tar.gz /data
```

### Restore Database

```bash
cat backup.sql | docker exec -i fidasim-postgres psql -U fidasim -d fidasim
```

## Updating

```bash
# Pull latest code
git pull

# Rebuild and restart
docker-compose -f docker-compose.prod.yaml build
docker-compose -f docker-compose.prod.yaml up -d

# Run database migrations if needed
docker exec fidasim-web-backend alembic upgrade head
```

## Documentation

- **DEPLOYMENT_GUIDE.md** - Detailed deployment instructions
- **QUICK_REFERENCE.md** - Command quick reference
- **INTEGRATION_FEATURE.md** - Visualization integration feature
- **API Docs** - http://localhost:8000/api/docs (when running)

## Support

For issues and questions:
1. Check logs: `docker-compose logs`
2. Review troubleshooting section above
3. Check GitHub issues
4. Review API documentation

## License

[Your License Here]

## Citation

If you use FIDASIM in your research, please cite:
[Citation information]

---

**Version**: 1.0.0
**Last Updated**: October 2025
**Status**: Production Ready ✅
