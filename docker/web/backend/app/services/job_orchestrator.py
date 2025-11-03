"""
Job Orchestrator Service

Manages communication between the web backend and fidasim-core API.
Handles job submission, status tracking, log retrieval, and cancellation.
"""

import httpx
from typing import Dict, Any, Optional, List
from uuid import UUID
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy.exc import SQLAlchemyError

from app.core.config import settings
from app.models.job import Job
from app.db.database import SessionLocal


class JobOrchestrator:
    """
    Orchestrates job execution by communicating with fidasim-core API
    """

    def __init__(self, db: Session):
        self.db = db
        self.core_api_url = settings.CORE_API_URL
        self.preprocessor_api_url = settings.PREPROCESSOR_API_URL

    async def create_job(
        self,
        user_id: UUID,
        run_id: str,
        config: Dict[str, Any],
        cores: int = 4,
        priority: int = 0
    ) -> Job:
        """
        Create a new job in the database and submit it to fidasim-core

        Args:
            user_id: User who submitted the job
            run_id: Unique run identifier
            config: Job configuration dictionary
            cores: Number of cores to allocate
            priority: Job priority (higher = more important)

        Returns:
            Created Job object
        """
        try:
            # Create job record
            job = Job(
                user_id=user_id,
                run_id=run_id,
                status="pending",
                config=config,
                cores_assigned=cores,
                priority=priority,
                input_files=config.get("input_files", []),
                output_files=[],
                error_message=None
            )

            self.db.add(job)
            self.db.commit()
            self.db.refresh(job)

            # Submit job to fidasim-core API in background (non-blocking)
            # Note: In production, this should be done via a task queue (Celery/Redis)
            # For now, we'll use asyncio to submit without blocking
            import asyncio
            asyncio.create_task(self._submit_to_core_async(job.id))

            return job

        except SQLAlchemyError as e:
            self.db.rollback()
            raise Exception(f"Database error creating job: {str(e)}")

    async def _submit_to_core_async(self, job_id: UUID):
        """
        Submit job to fidasim-core API asynchronously (non-blocking)

        Args:
            job_id: Job UUID
        """
        # Create new DB session for background task
        db = SessionLocal()
        try:
            job = db.query(Job).filter(Job.id == job_id).first()
            if not job:
                return

            await self._submit_to_core_impl(job, db)
        except Exception as e:
            import traceback
            error_detail = str(e) if str(e) else repr(e)
            if job:
                job.status = "failed"
                job.error_message = f"Failed to submit to core API: {error_detail}"
                db.commit()
            print(f"ERROR submitting job {job.run_id if job else job_id}: {error_detail}")
            print(f"Traceback: {traceback.format_exc()}")
        finally:
            db.close()

    async def _submit_to_core_impl(self, job: Job, db: Session):
        """
        Submit job to fidasim-core API

        Args:
            job: Job object to submit
        """
        async with httpx.AsyncClient(timeout=300.0) as client:  # 5 minute timeout for large file copies
            # Get input file path - check if it's a full path or just filename
            input_file = job.config.get("input_file", "")
            if input_file and not input_file.startswith("/"):
                # If just a filename, prepend the user's upload directory
                input_file = f"/data/uploads/{job.user_id}/{input_file}"

            payload = {
                "runid": job.run_id,
                "input_file": input_file,
                "cores": job.cores_assigned,
                "device": job.config.get("device", "fidasim_device"),
                "result_dir": job.config.get("result_dir", f"/data/results/{job.run_id}"),
                "tables_file": job.config.get("tables_file", "/app/tables/atomic_tables.h5"),
                "passive_only": job.config.get("passive_only", False)
            }

            response = await client.post(
                f"{self.core_api_url}/run",
                json=payload
            )

            if response.status_code == 200:
                result = response.json()
                # Update job with fidasim-core job_id
                job.metadata_ = job.metadata_ or {}
                job.metadata_["core_job_id"] = result.get("job_id")
                job.status = "running"
                job.started_at = datetime.utcnow()
                db.commit()  # Fixed: Use the correct db session
            else:
                raise Exception(f"Core API returned status {response.status_code}: {response.text}")

    async def get_job_status(self, job_id: UUID) -> Dict[str, Any]:
        """
        Get job status from fidasim-core API

        Args:
            job_id: Job UUID

        Returns:
            Job status information
        """
        job = self.db.query(Job).filter(Job.id == job_id).first()
        if not job:
            raise ValueError(f"Job {job_id} not found")

        # If job is already in terminal state, return cached status
        if job.status in ["completed", "failed", "cancelled"]:
            return {
                "job_id": str(job.id),
                "run_id": job.run_id,
                "status": job.status,
                "error_message": job.error_message,
                "output_files": job.output_files
            }

        # Get status from fidasim-core API
        core_job_id = job.metadata_.get("core_job_id") if job.metadata_ else None
        if not core_job_id:
            # Job not yet submitted to core
            return {
                "job_id": str(job.id),
                "run_id": job.run_id,
                "status": "pending",
                "error_message": None,
                "output_files": []
            }

        try:
            async with httpx.AsyncClient(timeout=30.0) as client:  # Increased: core API can be slow when FIDASIM is running
                response = await client.get(
                    f"{self.core_api_url}/status/{core_job_id}"
                )

                if response.status_code == 200:
                    core_status = response.json()

                    # Update job status from core API
                    job.status = core_status.get("status", job.status)

                    if job.status == "completed":
                        job.completed_at = datetime.utcnow()
                        job.output_files = core_status.get("output_files", [])
                    elif job.status == "failed":
                        job.error_message = core_status.get("error_message", "Unknown error")
                        job.completed_at = datetime.utcnow()

                    self.db.commit()

                    return core_status
                else:
                    # If we can't get status, return what we have
                    return {
                        "job_id": str(job.id),
                        "run_id": job.run_id,
                        "status": job.status,
                        "error_message": f"Failed to get status from core API: {response.status_code}",
                        "output_files": job.output_files
                    }

        except Exception as e:
            return {
                "job_id": str(job.id),
                "run_id": job.run_id,
                "status": job.status,
                "error_message": f"Error getting status: {str(e)}",
                "output_files": job.output_files
            }

    async def get_job_logs(self, job_id: UUID, tail: int = 100) -> Dict[str, Any]:
        """
        Get job logs from fidasim-core API

        Args:
            job_id: Job UUID
            tail: Number of lines to return

        Returns:
            Job logs
        """
        job = self.db.query(Job).filter(Job.id == job_id).first()
        if not job:
            raise ValueError(f"Job {job_id} not found")

        core_job_id = job.metadata_.get("core_job_id") if job.metadata_ else None
        if not core_job_id:
            return {
                "job_id": str(job.id),
                "run_id": job.run_id,
                "logs": "Job is being submitted to core API. Please wait a moment and refresh."
            }

        try:
            async with httpx.AsyncClient(timeout=30.0) as client:  # Increased: core API can be slow when FIDASIM is running
                response = await client.get(
                    f"{self.core_api_url}/logs/{core_job_id}",
                    params={"tail": tail}
                )

                if response.status_code == 200:
                    log_data = response.json()
                    # Core API returns {"lines": [...], "total_lines": N}
                    # Convert to plain text for display
                    if "lines" in log_data:
                        logs_text = ''.join(log_data["lines"])
                        return {
                            "job_id": str(job.id),
                            "run_id": job.run_id,
                            "logs": logs_text
                        }
                    else:
                        # Fallback for old format
                        return log_data
                else:
                    return {
                        "job_id": str(job.id),
                        "run_id": job.run_id,
                        "logs": f"Failed to retrieve logs: {response.status_code}"
                    }

        except Exception as e:
            return {
                "job_id": str(job.id),
                "run_id": job.run_id,
                "logs": f"Error retrieving logs: {str(e)}"
            }

    async def cancel_job(self, job_id: UUID) -> bool:
        """
        Cancel a running or pending job

        Args:
            job_id: Job UUID

        Returns:
            True if successfully cancelled
        """
        job = self.db.query(Job).filter(Job.id == job_id).first()
        if not job:
            raise ValueError(f"Job {job_id} not found")

        if job.status not in ["pending", "running"]:
            raise ValueError(f"Cannot cancel job with status '{job.status}'")

        core_job_id = job.metadata_.get("core_job_id") if job.metadata_ else None

        # If job was submitted to core, cancel it there
        if core_job_id:
            try:
                async with httpx.AsyncClient(timeout=10.0) as client:
                    response = await client.delete(
                        f"{self.core_api_url}/jobs/{core_job_id}/cancel"
                    )

                    if response.status_code not in [200, 204]:
                        raise Exception(f"Core API returned status {response.status_code}")

            except Exception as e:
                # Even if core cancellation fails, we'll mark it cancelled locally
                job.error_message = f"Warning: Core API cancellation failed: {str(e)}"

        # Update job status
        job.status = "cancelled"
        job.completed_at = datetime.utcnow()
        self.db.commit()

        return True

    async def get_job_results(self, job_id: UUID, filename: str) -> Optional[bytes]:
        """
        Download result file for a completed job

        Args:
            job_id: Job UUID
            filename: Name of the file to download

        Returns:
            File contents as bytes
        """
        job = self.db.query(Job).filter(Job.id == job_id).first()
        if not job:
            raise ValueError(f"Job {job_id} not found")

        if job.status != "completed":
            raise ValueError(f"Job is not completed (status: {job.status})")

        core_job_id = job.metadata_.get("core_job_id") if job.metadata_ else None
        if not core_job_id:
            raise ValueError("Job has no core job ID")

        try:
            async with httpx.AsyncClient(timeout=60.0) as client:
                response = await client.get(
                    f"{self.core_api_url}/results/{core_job_id}/{filename}"
                )

                if response.status_code == 200:
                    return response.content
                else:
                    raise Exception(f"Failed to download file: {response.status_code}")

        except Exception as e:
            raise Exception(f"Error downloading result file: {str(e)}")

    async def list_all_jobs(self, limit: int = 100, offset: int = 0) -> List[Job]:
        """
        List all jobs (admin function)

        Args:
            limit: Maximum number of jobs to return
            offset: Number of jobs to skip

        Returns:
            List of Job objects
        """
        jobs = self.db.query(Job).order_by(Job.created_at.desc()).offset(offset).limit(limit).all()
        return jobs

    async def sync_job_status(self, job_id: UUID):
        """
        Sync job status with fidasim-core API
        This can be called periodically to update job statuses

        Args:
            job_id: Job UUID
        """
        await self.get_job_status(job_id)
