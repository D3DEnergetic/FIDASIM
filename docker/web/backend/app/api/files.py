"""
Files API endpoints

Handles file upload, download, and management for FIDASIM jobs
"""

import os
import shutil
import aiofiles
from pathlib import Path
from typing import List, Optional
from uuid import UUID, uuid4
from datetime import datetime

from fastapi import (
    APIRouter,
    Depends,
    HTTPException,
    status,
    UploadFile,
    File,
    Query
)
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from app.core.security import get_current_active_user
from app.core.config import settings
from app.db.database import get_db
from app.models.user import User

router = APIRouter()

# Ensure upload directory exists
UPLOAD_DIR = Path(settings.UPLOAD_DIR)
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


def get_user_upload_dir(user_id: UUID) -> Path:
    """
    Get user-specific upload directory

    Args:
        user_id: User UUID

    Returns:
        Path to user's upload directory
    """
    user_dir = UPLOAD_DIR / str(user_id)
    user_dir.mkdir(parents=True, exist_ok=True)
    return user_dir


def sanitize_filename(filename: str) -> str:
    """
    Sanitize filename to prevent directory traversal attacks

    Args:
        filename: Original filename

    Returns:
        Sanitized filename
    """
    # Remove path components
    filename = os.path.basename(filename)
    # Remove potentially dangerous characters
    filename = filename.replace("..", "")
    return filename


@router.post("/upload")
async def upload_file(
    file: UploadFile = File(...),
    description: Optional[str] = None,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Upload a file for use in FIDASIM jobs

    Args:
        file: File to upload
        description: Optional file description
        current_user: Current authenticated user
        db: Database session

    Returns:
        File information
    """
    # Check file size
    file_content = await file.read()
    file_size = len(file_content)

    if file_size > settings.MAX_UPLOAD_SIZE:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"File too large. Maximum size is {settings.MAX_UPLOAD_SIZE / (1024**2):.0f} MB"
        )

    # Sanitize filename
    original_filename = sanitize_filename(file.filename)

    # Create unique filename to avoid collisions
    file_id = uuid4()
    file_extension = Path(original_filename).suffix
    unique_filename = f"{file_id}{file_extension}"

    # Get user directory
    user_dir = get_user_upload_dir(current_user.id)
    file_path = user_dir / unique_filename

    # Save file
    try:
        async with aiofiles.open(file_path, 'wb') as f:
            await f.write(file_content)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to save file: {str(e)}"
        )

    # Return file information
    return {
        "file_id": str(file_id),
        "filename": original_filename,
        "unique_filename": unique_filename,
        "size": file_size,
        "upload_date": datetime.utcnow().isoformat(),
        "path": str(file_path),
        "description": description
    }


@router.post("/upload-multiple")
async def upload_multiple_files(
    files: List[UploadFile] = File(...),
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Upload multiple files at once

    Args:
        files: List of files to upload
        current_user: Current authenticated user
        db: Database session

    Returns:
        List of uploaded file information
    """
    if len(files) > 20:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Maximum 20 files can be uploaded at once"
        )

    uploaded_files = []

    for file in files:
        try:
            result = await upload_file(file, None, current_user, db)
            uploaded_files.append(result)
        except HTTPException as e:
            # Return partial success with error information
            uploaded_files.append({
                "filename": file.filename,
                "error": e.detail
            })

    return {"files": uploaded_files}


@router.get("/list")
async def list_files(
    current_user: User = Depends(get_current_active_user)
):
    """
    List all uploaded files for current user

    Args:
        current_user: Current authenticated user

    Returns:
        List of files
    """
    user_dir = get_user_upload_dir(current_user.id)

    if not user_dir.exists():
        return {"files": []}

    files = []
    for file_path in user_dir.iterdir():
        if file_path.is_file():
            stat = file_path.stat()
            files.append({
                "filename": file_path.name,
                "size": stat.st_size,
                "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
                "path": str(file_path)
            })

    # Sort by modification time (newest first)
    files.sort(key=lambda x: x["modified"], reverse=True)

    return {"files": files}


@router.get("/download/{filename}")
async def download_file(
    filename: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Download a previously uploaded file

    Args:
        filename: Name of the file to download
        current_user: Current authenticated user

    Returns:
        File download
    """
    # Sanitize filename
    filename = sanitize_filename(filename)

    # Get user directory
    user_dir = get_user_upload_dir(current_user.id)
    file_path = user_dir / filename

    # Check if file exists
    if not file_path.exists():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    # Check if path is within user directory (security)
    try:
        file_path.resolve().relative_to(user_dir.resolve())
    except ValueError:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Access denied"
        )

    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type="application/octet-stream"
    )


@router.delete("/{filename}")
async def delete_file(
    filename: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Delete an uploaded file

    Args:
        filename: Name of the file to delete
        current_user: Current authenticated user

    Returns:
        Success message
    """
    # Sanitize filename
    filename = sanitize_filename(filename)

    # Get user directory
    user_dir = get_user_upload_dir(current_user.id)
    file_path = user_dir / filename

    # Check if file exists
    if not file_path.exists():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    # Check if path is within user directory (security)
    try:
        file_path.resolve().relative_to(user_dir.resolve())
    except ValueError:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Access denied"
        )

    # Delete file
    try:
        file_path.unlink()
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to delete file: {str(e)}"
        )

    return {"message": f"File '{filename}' deleted successfully"}


@router.get("/info/{filename}")
async def get_file_info(
    filename: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Get information about an uploaded file

    Args:
        filename: Name of the file
        current_user: Current authenticated user

    Returns:
        File information
    """
    # Sanitize filename
    filename = sanitize_filename(filename)

    # Get user directory
    user_dir = get_user_upload_dir(current_user.id)
    file_path = user_dir / filename

    # Check if file exists
    if not file_path.exists():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    # Check if path is within user directory (security)
    try:
        file_path.resolve().relative_to(user_dir.resolve())
    except ValueError:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Access denied"
        )

    # Get file info
    stat = file_path.stat()

    return {
        "filename": filename,
        "size": stat.st_size,
        "size_human": f"{stat.st_size / (1024**2):.2f} MB" if stat.st_size > 1024**2 else f"{stat.st_size / 1024:.2f} KB",
        "created": datetime.fromtimestamp(stat.st_ctime).isoformat(),
        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "path": str(file_path)
    }


@router.get("/storage-usage")
async def get_storage_usage(
    current_user: User = Depends(get_current_active_user)
):
    """
    Get storage usage statistics for current user

    Args:
        current_user: Current authenticated user

    Returns:
        Storage usage information
    """
    user_dir = get_user_upload_dir(current_user.id)

    if not user_dir.exists():
        return {
            "total_size": 0,
            "total_size_human": "0 B",
            "file_count": 0
        }

    total_size = 0
    file_count = 0

    for file_path in user_dir.rglob("*"):
        if file_path.is_file():
            total_size += file_path.stat().st_size
            file_count += 1

    # Format size
    if total_size > 1024**3:  # GB
        size_human = f"{total_size / (1024**3):.2f} GB"
    elif total_size > 1024**2:  # MB
        size_human = f"{total_size / (1024**2):.2f} MB"
    elif total_size > 1024:  # KB
        size_human = f"{total_size / 1024:.2f} KB"
    else:
        size_human = f"{total_size} B"

    return {
        "total_size": total_size,
        "total_size_human": size_human,
        "file_count": file_count,
        "max_size": settings.MAX_UPLOAD_SIZE,
        "max_size_human": f"{settings.MAX_UPLOAD_SIZE / (1024**3):.2f} GB"
    }
