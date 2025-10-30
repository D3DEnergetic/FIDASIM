"""
API endpoints for data visualization
"""

from fastapi import APIRouter, HTTPException, Depends, Query
from typing import Dict, List, Optional, Any
from pydantic import BaseModel
from sqlalchemy.orm import Session
from uuid import UUID
import logging

from ..core.security import get_current_active_user
from ..models.user import User
from ..models.job import Job
from ..db.database import get_db
from ..services.visualization import visualization_service

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/visualization", tags=["visualization"])

class SliceRequest(BaseModel):
    """Request model for data slicing"""
    variable: str
    slices: Optional[Dict[str, Any]] = None

class PlotRequest(BaseModel):
    """Request model for plot data"""
    variable: str
    x_dim: str
    y_dim: Optional[str] = None
    fixed_coords: Optional[Dict[str, int]] = None
    integrate_dims: Optional[List[str]] = None

@router.get("/files/{job_id}")
async def list_visualization_files(
    job_id: str,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
) -> List[Dict[str, Any]]:
    """
    List available output files for visualization

    Args:
        job_id: Job ID (UUID)

    Returns:
        List of available files
    """
    try:
        # Get the job from database to retrieve core_job_id
        job = db.query(Job).filter(Job.id == job_id).first()
        if not job:
            raise HTTPException(status_code=404, detail="Job not found")

        # Get core_job_id from metadata
        core_job_id = job.metadata_.get("core_job_id") if job.metadata_ else None
        if not core_job_id:
            # If no core_job_id, try using run_id as fallback
            core_job_id = job.run_id

        # Use core_job_id to list files
        files = visualization_service.list_available_files(core_job_id)
        return files
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error listing files for job {job_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/dataset/info")
async def get_dataset_info(
    file_path: str = Query(..., description="Path to the dataset file"),
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, Any]:
    """
    Get information about a dataset

    Args:
        file_path: Path to the dataset file

    Returns:
        Dataset metadata including dimensions, variables, and attributes
    """
    try:
        info = visualization_service.get_dataset_info(file_path)
        return info
    except Exception as e:
        logger.error(f"Error getting dataset info for {file_path}: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/dataset/variable")
async def get_variable_data(
    file_path: str,
    request: SliceRequest,
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, Any]:
    """
    Get data for a specific variable with optional slicing

    Args:
        file_path: Path to the dataset
        request: Variable name and optional slices

    Returns:
        Variable data and metadata
    """
    try:
        data = visualization_service.get_variable_data(
            file_path=file_path,
            variable=request.variable,
            slices=request.slices
        )
        return data
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting variable data: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/plot/data")
async def get_plot_data(
    file_path: str,
    request: PlotRequest,
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, Any]:
    """
    Get data formatted for plotting

    Args:
        file_path: Path to the dataset
        request: Plot configuration

    Returns:
        Plot-ready data
    """
    try:
        data = visualization_service.get_slice_for_plot(
            file_path=file_path,
            variable=request.variable,
            x_dim=request.x_dim,
            y_dim=request.y_dim,
            fixed_coords=request.fixed_coords,
            integrate_dims=request.integrate_dims
        )
        return data
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting plot data: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/dataset/variables")
async def list_variables(
    file_path: str = Query(..., description="Path to the dataset file"),
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, List[str]]:
    """
    List all variables and coordinates in a dataset

    Args:
        file_path: Path to the dataset

    Returns:
        Lists of data variables and coordinates
    """
    try:
        info = visualization_service.get_dataset_info(file_path)
        return {
            'data_variables': list(info['data_variables'].keys()),
            'coordinates': list(info['coordinates'].keys()),
            'dimensions': list(info['dimensions'].keys())
        }
    except Exception as e:
        logger.error(f"Error listing variables: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/cache/clear")
async def clear_cache(
    file_path: Optional[str] = None,
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, str]:
    """
    Clear cached datasets

    Args:
        file_path: Optional specific file to clear from cache

    Returns:
        Success message
    """
    try:
        visualization_service.clear_cache(file_path)
        message = f"Cache cleared for {file_path}" if file_path else "All cache cleared"
        return {"message": message}
    except Exception as e:
        logger.error(f"Error clearing cache: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))