"""
Job Status Synchronization Service

Periodically syncs job statuses from fidasim-core API to web backend database
"""

import asyncio
import logging
from datetime import datetime
from sqlalchemy.orm import Session

from app.db.database import SessionLocal
from app.models.job import Job
from app.services.job_orchestrator import JobOrchestrator

logger = logging.getLogger(__name__)


async def sync_running_jobs():
    """
    Sync status for all jobs marked as 'running' or 'pending'
    """
    db: Session = SessionLocal()

    try:
        # Get all running/pending jobs
        running_jobs = db.query(Job).filter(
            Job.status.in_(["running", "pending"])
        ).all()

        if not running_jobs:
            logger.debug("No running/pending jobs to sync")
            return

        logger.info(f"Syncing status for {len(running_jobs)} running/pending jobs")

        orchestrator = JobOrchestrator(db)

        for job in running_jobs:
            try:
                # Sync status from core API
                await orchestrator.sync_job_status(job.id)
                logger.info(f"Synced status for job {job.run_id}")
            except Exception as e:
                logger.error(f"Failed to sync job {job.run_id}: {e}")

    except Exception as e:
        logger.error(f"Error in job sync: {e}")
    finally:
        db.close()


async def job_sync_worker(interval: int = 10):
    """
    Background worker that syncs job statuses periodically

    Args:
        interval: Sync interval in seconds (default: 10)
    """
    logger.info(f"Job sync worker started (interval: {interval}s)")

    while True:
        try:
            await sync_running_jobs()
        except Exception as e:
            logger.error(f"Error in job sync worker: {e}")

        # Wait before next sync
        await asyncio.sleep(interval)


def start_job_sync_worker(interval: int = 10):
    """
    Start the job sync worker in the background

    Args:
        interval: Sync interval in seconds (default: 10)
    """
    asyncio.create_task(job_sync_worker(interval))
    logger.info("Job sync worker task created")
