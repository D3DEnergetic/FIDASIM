"""
FIDASIM Core API Service
Provides REST endpoints for running FIDASIM simulations
"""

import os
import asyncio
import subprocess
import uuid
import json
import logging
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any
from enum import Enum

from fastapi import FastAPI, HTTPException, BackgroundTasks, UploadFile, File, Form
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel, Field
import h5py
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = Path(os.getenv("DATA_DIR", "/data"))
FIDASIM_EXECUTABLE = Path("/app/fidasim")
ATOMIC_TABLES = Path("/app/tables/atomic_tables.h5")
MAX_CORES = int(os.getenv("MAX_CORES", "32"))

# Ensure data directory exists
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Job tracking (in production, this would be in Redis/Database)
jobs: Dict[str, Dict] = {}

class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

class SimulationConfig(BaseModel):
    """Configuration for a FIDASIM simulation run"""
    runid: str = Field(description="Unique run identifier")
    input_file: str = Field(description="Path to input namelist file")
    cores: int = Field(default=4, ge=1, le=32, description="Number of OpenMP threads")
    equilibrium_file: Optional[str] = Field(default=None, description="Path to equilibrium HDF5 file")
    geometry_file: Optional[str] = Field(default=None, description="Path to geometry HDF5 file")
    distribution_file: Optional[str] = Field(default=None, description="Path to distribution HDF5 file")
    neutrals_file: Optional[str] = Field(default=None, description="Path to neutrals HDF5 file")

class JobInfo(BaseModel):
    """Information about a simulation job"""
    job_id: str
    run_id: str
    status: JobStatus
    cores: int
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    output_files: List[str] = []
    log_file: Optional[str] = None

class ValidationResult(BaseModel):
    """Result of input file validation"""
    valid: bool
    errors: List[str] = []
    warnings: List[str] = []
    info: Dict[str, Any] = {}

# Initialize FastAPI app
app = FastAPI(
    title="FIDASIM Core API",
    description="API for running FIDASIM simulations",
    version="1.0.0"
)

@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "service": "FIDASIM Core",
        "version": "3.0.0-dev",
        "status": "operational",
        "executable": str(FIDASIM_EXECUTABLE),
        "max_cores": MAX_CORES
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    checks = {
        "api": True,
        "executable_exists": FIDASIM_EXECUTABLE.exists(),
        "atomic_tables_exists": ATOMIC_TABLES.exists(),
        "data_dir_writable": os.access(DATA_DIR, os.W_OK)
    }

    healthy = all(checks.values())

    return {
        "healthy": healthy,
        "checks": checks,
        "timestamp": datetime.utcnow().isoformat()
    }

@app.post("/validate")
async def validate_inputs(
    input_file: UploadFile = File(...),
    equilibrium_file: Optional[UploadFile] = File(None),
    geometry_file: Optional[UploadFile] = File(None),
    distribution_file: Optional[UploadFile] = File(None)
) -> ValidationResult:
    """Validate input files before running simulation"""

    errors = []
    warnings = []
    info = {}

    # Create temporary directory for validation
    temp_dir = DATA_DIR / f"validate_{uuid.uuid4().hex}"
    temp_dir.mkdir(parents=True)

    try:
        # Save input namelist
        input_path = temp_dir / "inputs.dat"
        content = await input_file.read()
        input_path.write_bytes(content)

        # Parse namelist to extract key parameters
        namelist_content = content.decode('utf-8')

        # Basic namelist validation
        if "&fidasim_inputs" not in namelist_content:
            errors.append("Missing &fidasim_inputs namelist group")

        # Check for required parameters
        required_params = ["shot", "time", "runid", "result_dir", "tables_file"]
        for param in required_params:
            if param not in namelist_content:
                warnings.append(f"Parameter '{param}' not found in namelist")

        # Extract simulation settings
        if "calc_bes" in namelist_content:
            info["calc_bes"] = "T" in namelist_content.split("calc_bes")[1][:20]
        if "calc_fida" in namelist_content:
            info["calc_fida"] = "T" in namelist_content.split("calc_fida")[1][:20]
        if "calc_npa" in namelist_content:
            info["calc_npa"] = "T" in namelist_content.split("calc_npa")[1][:20]

        # Validate HDF5 files if provided
        hdf5_files = {
            "equilibrium": equilibrium_file,
            "geometry": geometry_file,
            "distribution": distribution_file
        }

        for file_type, file_obj in hdf5_files.items():
            if file_obj:
                file_path = temp_dir / f"{file_type}.h5"
                content = await file_obj.read()
                file_path.write_bytes(content)

                try:
                    with h5py.File(file_path, 'r') as f:
                        info[f"{file_type}_groups"] = list(f.keys())[:10]  # First 10 groups
                except Exception as e:
                    errors.append(f"Invalid {file_type} HDF5 file: {str(e)}")

        valid = len(errors) == 0

    except Exception as e:
        errors.append(f"Validation error: {str(e)}")
        valid = False

    finally:
        # Clean up temporary files
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    return ValidationResult(
        valid=valid,
        errors=errors,
        warnings=warnings,
        info=info
    )

async def run_fidasim_task(job_id: str, config: SimulationConfig):
    """Background task to run FIDASIM simulation"""

    job = jobs[job_id]
    job_dir = DATA_DIR / f"job_{job_id}"

    try:
        # Update job status
        job["status"] = JobStatus.RUNNING
        job["started_at"] = datetime.utcnow()

        # Set up environment - ensure OMP_NUM_THREADS is set correctly
        env = {}
        env["OMP_NUM_THREADS"] = str(config.cores)
        env["PATH"] = os.environ.get("PATH", "/usr/local/bin:/usr/bin:/bin")
        env["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "/usr/local/lib")

        logger.info(f"Setting OMP_NUM_THREADS={config.cores} for job {job_id}")

        # Prepare input file path
        input_path = job_dir / Path(config.input_file).name

        # Rewrite tables_file path in the input file to use container's tables
        # Read the input file and replace tables_file path
        with open(input_path, 'r') as f:
            input_content = f.read()

        # Replace file paths to use container paths
        import re
        # Replace tables_file with container's atomic tables
        input_content = re.sub(
            r"(tables_file\s*=\s*)'[^']*'",
            r"\1'/app/tables/atomic_tables.h5'",
            input_content
        )
        # Replace result_dir to use job directory so outputs go there
        input_content = re.sub(
            r"(result_dir\s*=\s*)'[^']*'",
            f"\\1'{job_dir}'",
            input_content
        )
        # Replace any test file paths to use /test mount point
        input_content = re.sub(
            r"'/home/prechelg/Work/FIDASIM/test/",
            "'/test/",
            input_content
        )
        input_content = re.sub(
            r"'/home/prechelg/FIDASIM/test/",
            "'/test/",
            input_content
        )

        # Copy required HDF5 files from /test to job directory if they exist
        # This avoids HDF5 permission issues with mounted read-only files
        test_files_to_copy = []
        for line in input_content.split('\n'):
            if '.h5' in line and 'file' in line and '=' in line:
                # Extract the file path from lines like: geometry_file = '/test/file.h5'
                match = re.search(r"=\s*'([^']*\.h5)'", line)
                if match:
                    filepath = match.group(1)
                    if filepath.startswith('/test/'):
                        filename = Path(filepath).name
                        src = Path(filepath)
                        dst = job_dir / filename
                        if src.exists():
                            shutil.copy2(src, dst)
                            # Update path in input content to use local copy
                            input_content = input_content.replace(filepath, str(dst))
                            logger.info(f"Copied {filename} to job directory")

        # Write modified input file
        with open(input_path, 'w') as f:
            f.write(input_content)

        # Create log file
        log_file = job_dir / f"{config.runid}.log"
        job["log_file"] = str(log_file)

        # Run FIDASIM
        logger.info(f"Starting FIDASIM for job {job_id} with {config.cores} cores (using /app/tables/atomic_tables.h5)")

        with open(log_file, 'w') as log:
            # FIDASIM command: /app/fidasim input_file.dat num_threads
            process = subprocess.Popen(
                [str(FIDASIM_EXECUTABLE), str(input_path), str(config.cores)],
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=str(job_dir)
            )

            # Store process ID for potential cancellation
            job["pid"] = process.pid

            # Wait for completion
            return_code = process.wait()

        if return_code == 0:
            job["status"] = JobStatus.COMPLETED
            logger.info(f"Job {job_id} completed successfully")

            # Find output files
            output_files = []
            for pattern in ["*_spectra.h5", "*_npa.h5", "*_neutrons.h5",
                           "*_birth.h5", "*_weights.h5"]:
                output_files.extend(job_dir.glob(pattern))

            job["output_files"] = [str(f.relative_to(job_dir)) for f in output_files]
        else:
            job["status"] = JobStatus.FAILED
            job["error_message"] = f"FIDASIM exited with code {return_code}"
            logger.error(f"Job {job_id} failed with return code {return_code}")

    except Exception as e:
        job["status"] = JobStatus.FAILED
        job["error_message"] = str(e)
        logger.error(f"Job {job_id} failed with error: {e}")

    finally:
        job["completed_at"] = datetime.utcnow()
        job.pop("pid", None)

@app.post("/run")
async def run_simulation(
    background_tasks: BackgroundTasks,
    config: SimulationConfig
) -> JobInfo:
    """Start a new FIDASIM simulation"""

    # Generate job ID
    job_id = uuid.uuid4().hex

    # Create job directory
    job_dir = DATA_DIR / f"job_{job_id}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # Validate input file exists
    input_path = DATA_DIR / config.input_file
    if not input_path.exists():
        raise HTTPException(status_code=404, detail=f"Input file not found: {config.input_file}")

    # Copy input file to job directory
    shutil.copy2(input_path, job_dir / input_path.name)

    # Copy other required files if specified
    for file_attr in ["equilibrium_file", "geometry_file", "distribution_file", "neutrals_file"]:
        file_path = getattr(config, file_attr)
        if file_path:
            src = DATA_DIR / file_path
            if src.exists():
                shutil.copy2(src, job_dir / src.name)
            else:
                logger.warning(f"File not found: {file_path}")

    # Create job entry
    job = {
        "job_id": job_id,
        "run_id": config.runid,
        "status": JobStatus.PENDING,
        "cores": config.cores,
        "created_at": datetime.utcnow(),
        "started_at": None,
        "completed_at": None,
        "error_message": None,
        "output_files": [],
        "log_file": None,
        "config": config.dict()
    }

    jobs[job_id] = job

    # Start background task
    background_tasks.add_task(run_fidasim_task, job_id, config)

    return JobInfo(**job)

@app.get("/jobs")
async def list_jobs(
    status: Optional[JobStatus] = None,
    limit: int = 100
) -> List[JobInfo]:
    """List all jobs, optionally filtered by status"""

    filtered_jobs = []
    for job in jobs.values():
        if status is None or job["status"] == status:
            filtered_jobs.append(JobInfo(**job))

    # Sort by creation time, most recent first
    filtered_jobs.sort(key=lambda x: x.created_at, reverse=True)

    return filtered_jobs[:limit]

@app.get("/status/{job_id}")
async def get_job_status(job_id: str) -> JobInfo:
    """Get the status of a specific job"""

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    return JobInfo(**jobs[job_id])

@app.get("/logs/{job_id}")
async def get_job_logs(
    job_id: str,
    tail: int = 100,
    stream: bool = False
):
    """Get logs for a specific job"""

    # Try to get log file from job data in memory
    log_file = None
    if job_id in jobs:
        log_file = jobs[job_id].get("log_file")

    # If not in memory, try to find the log file on disk
    if not log_file:
        job_dir = DATA_DIR / f"job_{job_id}"
        if job_dir.exists():
            # Find any .log file in the job directory
            log_files = list(job_dir.glob("*.log"))
            if log_files:
                log_file = str(log_files[0])

    if not log_file or not Path(log_file).exists():
        return {"lines": [], "total_lines": 0}

    if stream and job["status"] == JobStatus.RUNNING:
        # Stream logs for running job
        async def log_streamer():
            with open(log_file, 'r') as f:
                # Go to end of file
                f.seek(0, 2)
                while job["status"] == JobStatus.RUNNING:
                    line = f.readline()
                    if line:
                        yield line
                    else:
                        await asyncio.sleep(0.5)
                # Send remaining lines
                for line in f:
                    yield line

        return StreamingResponse(log_streamer(), media_type="text/plain")
    else:
        # Return tail of log file
        with open(log_file, 'r') as f:
            lines = f.readlines()
            return {"lines": lines[-tail:], "total_lines": len(lines)}

@app.delete("/jobs/{job_id}")
async def cancel_job(job_id: str):
    """Cancel a running job"""

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    job = jobs[job_id]

    if job["status"] != JobStatus.RUNNING:
        raise HTTPException(
            status_code=400,
            detail=f"Job is not running (status: {job['status']})"
        )

    # Kill the process if it's running
    pid = job.get("pid")
    if pid:
        try:
            os.kill(pid, 15)  # SIGTERM
            job["status"] = JobStatus.CANCELLED
            job["completed_at"] = datetime.utcnow()
            return {"message": f"Job {job_id} cancelled"}
        except ProcessLookupError:
            job["status"] = JobStatus.FAILED
            job["error_message"] = "Process not found"

    raise HTTPException(status_code=500, detail="Failed to cancel job")

@app.get("/results/{job_id}/{filename}")
async def download_result(job_id: str, filename: str):
    """Download a result file from a completed job"""

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")

    job = jobs[job_id]

    if job["status"] not in [JobStatus.COMPLETED, JobStatus.FAILED]:
        raise HTTPException(
            status_code=400,
            detail=f"Job not completed (status: {job['status']})"
        )

    # Security check: ensure filename is in output_files
    if filename not in job["output_files"] and filename != f"{job['run_id']}.log":
        raise HTTPException(status_code=404, detail="File not found")

    file_path = DATA_DIR / f"job_{job_id}" / filename

    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type="application/octet-stream"
    )

@app.post("/upload")
async def upload_file(
    file: UploadFile = File(...),
    directory: str = Form("")
):
    """Upload a file to the data directory"""

    # Sanitize directory name
    if directory:
        upload_dir = DATA_DIR / directory
        upload_dir.mkdir(parents=True, exist_ok=True)
    else:
        upload_dir = DATA_DIR

    # Sanitize filename
    filename = Path(file.filename).name
    file_path = upload_dir / filename

    # Save file
    content = await file.read()
    file_path.write_bytes(content)

    return {
        "filename": filename,
        "path": str(file_path.relative_to(DATA_DIR)),
        "size": len(content)
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)