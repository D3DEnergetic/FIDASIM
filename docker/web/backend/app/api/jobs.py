"""
Jobs API endpoints
"""

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, and_
from typing import List, Optional
from uuid import UUID

from app.core.security import get_current_active_user
from app.db.database import get_db
from app.models.user import User
from app.models.job import Job
from app.schemas.job import JobSubmit, JobResponse, JobDetail, JobStats, JobUpdate
from app.services.job_orchestrator import JobOrchestrator

router = APIRouter()


@router.post("/submit", response_model=JobResponse, status_code=status.HTTP_201_CREATED)
async def submit_job(
    job_data: JobSubmit,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Submit a new FIDASIM job
    """
    # Check if run_id already exists
    existing_job = db.query(Job).filter(Job.run_id == job_data.run_id).first()
    if existing_job:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Job with run_id '{job_data.run_id}' already exists"
        )

    # Check user's job limit
    active_jobs = db.query(Job).filter(
        and_(
            Job.user_id == current_user.id,
            Job.status.in_(["pending", "running"])
        )
    ).count()

    if active_jobs >= current_user.max_jobs:
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=f"Maximum concurrent jobs limit ({current_user.max_jobs}) reached"
        )

    # Create job
    orchestrator = JobOrchestrator(db)
    job = await orchestrator.create_job(
        user_id=current_user.id,
        run_id=job_data.run_id,
        config=job_data.config.dict(),
        cores=job_data.config.cores,
        priority=job_data.priority
    )

    return job


@router.get("/", response_model=List[JobResponse])
async def list_jobs(
    status_filter: Optional[str] = Query(None, regex="^(pending|running|completed|failed|cancelled)$"),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    List jobs for current user
    """
    query = db.query(Job).filter(Job.user_id == current_user.id)

    if status_filter:
        query = query.filter(Job.status == status_filter)

    query = query.order_by(Job.created_at.desc())
    jobs = query.offset(offset).limit(limit).all()

    return jobs


@router.get("/stats", response_model=JobStats)
async def get_job_stats(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Get job statistics for current user
    """
    user_jobs = db.query(Job).filter(Job.user_id == current_user.id)

    stats = {
        "total": user_jobs.count(),
        "pending": user_jobs.filter(Job.status == "pending").count(),
        "running": user_jobs.filter(Job.status == "running").count(),
        "completed": user_jobs.filter(Job.status == "completed").count(),
        "failed": user_jobs.filter(Job.status == "failed").count(),
        "cancelled": user_jobs.filter(Job.status == "cancelled").count(),
    }

    return stats


@router.get("/{job_id}", response_model=JobDetail)
async def get_job(
    job_id: UUID,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Get detailed information about a specific job
    """
    job = db.query(Job).filter(Job.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )

    # Check ownership
    if job.user_id != current_user.id and current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to access this job"
        )

    return job


@router.patch("/{job_id}", response_model=JobResponse)
async def update_job(
    job_id: UUID,
    job_update: JobUpdate,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Update a job (admin or owner only)
    """
    job = db.query(Job).filter(Job.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )

    # Check ownership
    if job.user_id != current_user.id and current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to modify this job"
        )

    # Update fields
    for field, value in job_update.dict(exclude_unset=True).items():
        setattr(job, field, value)

    db.commit()
    db.refresh(job)

    return job


@router.delete("/{job_id}/cancel")
async def cancel_job(
    job_id: UUID,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Cancel a running or pending job
    """
    job = db.query(Job).filter(Job.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )

    # Check ownership
    if job.user_id != current_user.id and current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to cancel this job"
        )

    if job.status not in ["pending", "running"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot cancel job with status '{job.status}'"
        )

    # Cancel the job
    orchestrator = JobOrchestrator(db)
    await orchestrator.cancel_job(job_id)

    return {"message": f"Job {job_id} cancelled successfully"}


@router.get("/{job_id}/logs")
async def get_job_logs(
    job_id: UUID,
    tail: int = Query(100, ge=1, le=10000),
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Get logs for a specific job
    """
    job = db.query(Job).filter(Job.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found"
        )

    # Check ownership
    if job.user_id != current_user.id and current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to access this job"
        )

    # Get logs from fidasim-core API
    orchestrator = JobOrchestrator(db)
    logs = await orchestrator.get_job_logs(job_id, tail)

    return logs