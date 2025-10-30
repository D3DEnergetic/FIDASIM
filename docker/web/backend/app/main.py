"""
FIDASIM Web API - Main Application
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from app.core.config import settings
from app.api import auth, jobs, users, files, websocket, visualization
from app.db.database import engine, Base


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Lifecycle events"""
    # Startup
    print("Starting FIDASIM Web API...")
    # Create database tables
    Base.metadata.create_all(bind=engine)
    print("Database tables created/verified")

    # Start background job sync worker
    from app.services.job_sync import start_job_sync_worker
    start_job_sync_worker(interval=10)  # Sync every 10 seconds
    print("Job sync worker started")

    yield

    # Shutdown
    print("Shutting down FIDASIM Web API...")


# Create FastAPI application
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="Web interface for FIDASIM - Fast-Ion Diagnostic Simulation",
    lifespan=lifespan,
    docs_url="/api/docs",
    redoc_url="/api/redoc",
    openapi_url="/api/openapi.json"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(auth.router, prefix="/api/auth", tags=["Authentication"])
app.include_router(users.router, prefix="/api/users", tags=["Users"])
app.include_router(jobs.router, prefix="/api/jobs", tags=["Jobs"])
app.include_router(files.router, prefix="/api/files", tags=["Files"])
app.include_router(websocket.router, prefix="/api/ws", tags=["WebSocket"])
app.include_router(visualization.router, prefix="/api", tags=["Visualization"])


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "operational",
        "docs": "/api/docs"
    }


@app.get("/api/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "fidasim-web",
        "version": settings.APP_VERSION
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="0.0.0.0", port=8000, reload=True)