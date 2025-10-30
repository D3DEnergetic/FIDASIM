# FIDASIM Web Platform - Essential Files for GitHub

This document lists the files and directories needed for someone to clone, build, and deploy the FIDASIM web platform.

## Essential Directories

### 1. `docker/` - Container Definitions
**Required:**
- `docker/web/backend/` - Complete backend application
  - `docker/web/backend/Dockerfile`
  - `docker/web/backend/requirements.txt`
  - `docker/web/backend/app/` - All Python application code
    - `docker/web/backend/app/main.py`
    - `docker/web/backend/app/api/` - API endpoints
    - `docker/web/backend/app/core/` - Core functionality (config, security)
    - `docker/web/backend/app/db/` - Database configuration
    - `docker/web/backend/app/models/` - Database models
    - `docker/web/backend/app/schemas/` - Pydantic schemas
    - `docker/web/backend/app/services/` - Business logic (including visualization)

- `docker/web/frontend/` - Complete frontend application
  - `docker/web/frontend/Dockerfile`
  - `docker/web/frontend/package.json`
  - `docker/web/frontend/package-lock.json` (if using npm ci)
  - `docker/web/frontend/tsconfig.json`
  - `docker/web/frontend/vite.config.ts`
  - `docker/web/frontend/index.html`
  - `docker/web/frontend/nginx.conf`
  - `docker/web/frontend/docker-entrypoint.sh`
  - `docker/web/frontend/src/` - All React application code
    - `docker/web/frontend/src/main.tsx`
    - `docker/web/frontend/src/App.tsx`
    - `docker/web/frontend/src/components/` - All React components
    - `docker/web/frontend/src/pages/` - All page components
    - `docker/web/frontend/src/services/` - API services
    - `docker/web/frontend/src/contexts/` - React contexts
    - `docker/web/frontend/src/types/` - TypeScript types
    - `docker/web/frontend/src/hooks/` - Custom hooks
    - `docker/web/frontend/src/utils/` - Utility functions

- `docker/core/` - FIDASIM core simulation service
  - `docker/core/Dockerfile`
  - `docker/core/api.py`
  - `docker/core/entrypoint.sh`
  - `docker/core/requirements.txt`

- `docker/preprocessor/` - Preprocessing service
  - `docker/preprocessor/Dockerfile`
  - `docker/preprocessor/api.py`
  - `docker/preprocessor/requirements.txt`

- `docker/postgres/` - Database initialization
  - `docker/postgres/init.sql` (if exists)

### 2. `scripts/` - Build and Deployment Scripts
**Required:**
- `scripts/build-web.sh` - Build web services
- `scripts/deploy-web.sh` - Deploy complete stack
- `scripts/test-web.sh` - Test deployment

**Optional but useful:**
- `scripts/build-all.sh` - Build all services
- `scripts/docker-diagnostic.sh` - Troubleshooting

### 3. Docker Compose Files
**Required:**
- `docker-compose.prod.yaml` - Production deployment configuration

**Optional:**
- `docker-compose.yaml` - Core services only
- `docker-compose.dev.yaml` - Development configuration

### 4. Environment Configuration
**Required:**
- `.env.prod.example` - Template for production environment variables

**Note:** Do NOT commit `.env` (actual secrets)

### 5. Documentation
**Highly Recommended:**
- `DEPLOYMENT_GUIDE.md` - Complete deployment instructions
- `QUICK_REFERENCE.md` - Quick command reference
- `README_WEB_PLATFORM.md` - Platform overview
- `DEPLOYMENT_SUCCESS.md` - Successful deployment notes

**Optional but useful:**
- `VISUALIZATION_CRASH_FIX.md` - Bug fix documentation
- `INTEGRATION_FEATURE.md` - Integration feature documentation
- `PROJECT_COMPLETE.md` - Project summary
- `FRONTEND_COMPLETE.md` - Frontend documentation
- `BACKEND_COMPLETE.md` - Backend documentation

## Files to Exclude (.gitignore)

Create or update `.gitignore` to exclude:

```gitignore
# Environment files with secrets
.env
.env.local

# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
*.egg-info/
dist/
build/

# Node
node_modules/
npm-debug.log*
yarn-debug.log*
yarn-error.log*
dist/
.vite/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db

# Docker volumes and data
/data/
*.log

# Test files (unless needed for testing)
test/
tests/

# Build artifacts
*.o
*.mod
*.a
fidasim

# Fortran source (if web-only deployment)
src/*.f90
deps/

# Temporary documentation
*_COMPLETE.md
*_STATUS.md
*_PROGRESS.md
*_SUMMARY.md
SESSION_*.md
CURRENT_STATUS.md
automation.md
```

## Minimal Deployment Package

For a minimal web-only deployment, you need:

```
FIDASIM/
├── docker/
│   ├── core/
│   ├── preprocessor/
│   ├── postgres/
│   └── web/
│       ├── backend/
│       └── frontend/
├── scripts/
│   ├── build-web.sh
│   ├── deploy-web.sh
│   └── test-web.sh
├── docker-compose.prod.yaml
├── .env.prod.example
├── DEPLOYMENT_GUIDE.md
├── QUICK_REFERENCE.md
├── README.md (or README_WEB_PLATFORM.md)
└── .gitignore
```

## NOT Needed for Web Platform Deployment

These are only needed if building FIDASIM from source:
- `src/` - Fortran source code
- `deps/` - Build dependencies
- `test/*.h5` - Test simulation files (though a few samples might be useful)
- Build scripts for Fortran compilation
- `.travis.yml` - CI configuration for core FIDASIM

## Repository Structure Options

### Option 1: Full Repository
Keep everything together, useful if you want to build FIDASIM from source and deploy the web platform.

### Option 2: Web Platform Branch
Create a separate branch with only web platform files:
```bash
git checkout -b web-platform
git rm -r src/ deps/ test/
git commit -m "Web platform deployment branch"
```

### Option 3: Separate Repository
Create a new repository with just the web platform components, linking to pre-built FIDASIM Docker images.

## Quick Start for New Users

After cloning, users should:

1. **Copy environment template:**
   ```bash
   cp .env.prod.example .env
   ```

2. **Configure environment:**
   ```bash
   # Edit .env and set:
   # - SECRET_KEY (generate with: openssl rand -hex 32)
   # - POSTGRES_PASSWORD
   # - VITE_API_URL (if remote deployment)
   # - VITE_WS_URL (if remote deployment)
   ```

3. **Deploy:**
   ```bash
   ./scripts/deploy-web.sh
   ```

4. **Create admin user:**
   ```bash
   docker exec fidasim-web-backend python -c "..."
   # (Full command in DEPLOYMENT_GUIDE.md)
   ```

5. **Access:**
   - Navigate to http://localhost:3000
   - Login with created credentials

## Docker Images

If you want to pre-build and publish Docker images:

```bash
# Build and tag
docker-compose -f docker-compose.prod.yaml build
docker tag fidasim-web-backend:latest yourregistry/fidasim-web-backend:v1.0
docker tag fidasim-web-frontend:latest yourregistry/fidasim-web-frontend:v1.0

# Push
docker push yourregistry/fidasim-web-backend:v1.0
docker push yourregistry/fidasim-web-frontend:v1.0
```

Then users can pull pre-built images instead of building from source.

## Summary

**Minimum for web platform deployment:**
- `docker/` directory (all subdirectories)
- `scripts/` directory (web-related scripts)
- `docker-compose.prod.yaml`
- `.env.prod.example`
- `DEPLOYMENT_GUIDE.md`
- `.gitignore`

**Total size estimate:** ~5-10 MB (excluding node_modules, which are downloaded during build)

**Build time:** 5-10 minutes on first build (with good internet connection)

**Deployment time:** 1-2 minutes after images are built
