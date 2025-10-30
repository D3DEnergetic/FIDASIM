"""
Job schemas for API validation
"""

from pydantic import BaseModel, Field
from typing import Optional, Dict, List, Any
from datetime import datetime
from uuid import UUID


class JobConfig(BaseModel):
    """Configuration for a FIDASIM job"""
    input_file: str = Field(..., description="Path to input file")
    device: Optional[str] = Field(default="fidasim_device", description="Device name")
    result_dir: Optional[str] = Field(None, description="Result directory")
    tables_file: Optional[str] = Field(default="/app/tables/atomic_tables.h5", description="Atomic tables file")
    passive_only: bool = Field(default=False, description="Run passive diagnostics only")
    cores: int = Field(default=4, ge=1, le=32, description="Number of cores")
    # Legacy fields
    equilibrium_file: Optional[str] = None
    geometry_file: Optional[str] = None
    distribution_file: Optional[str] = None
    neutrals_file: Optional[str] = None
    additional_params: Optional[Dict[str, Any]] = Field(default_factory=dict)


class JobSubmit(BaseModel):
    """Schema for submitting a new job"""
    run_id: str = Field(..., min_length=1, max_length=255)
    config: JobConfig
    priority: int = Field(default=0, ge=0, le=10)


class JobUpdate(BaseModel):
    """Schema for updating a job"""
    status: Optional[str] = Field(None, pattern="^(pending|running|completed|failed|cancelled)$")
    priority: Optional[int] = Field(None, ge=1, le=10)


class JobResponse(BaseModel):
    """Job response schema"""
    id: UUID
    user_id: Optional[UUID]
    run_id: str
    status: str
    cores_assigned: int
    priority: int
    created_at: datetime
    queued_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    output_files: List[str] = []
    error_message: Optional[str] = None

    class Config:
        from_attributes = True


class JobDetail(JobResponse):
    """Detailed job information"""
    input_files: List[str] = []  # Changed from Dict to List to match model
    config: Dict[str, Any]
    container_id: Optional[str] = None
    exit_code: Optional[int] = None
    cpu_time_seconds: Optional[float] = None
    memory_peak_mb: Optional[float] = None
    particles_simulated: Optional[int] = None


class JobStats(BaseModel):
    """Job statistics"""
    total: int
    pending: int
    running: int
    completed: int
    failed: int
    cancelled: int