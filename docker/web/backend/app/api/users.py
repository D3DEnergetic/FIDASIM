"""
Users API endpoints

Provides user management functionality (admin only)
"""

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from uuid import UUID
from datetime import datetime

from app.core.security import (
    get_current_active_user,
    require_admin,
    get_password_hash
)
from app.db.database import get_db
from app.models.user import User
from app.schemas.user import UserCreate, UserUpdate, UserResponse

router = APIRouter()


@router.get("/", response_model=List[UserResponse])
async def list_users(
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    role_filter: Optional[str] = Query(None, regex="^(admin|user|viewer)$"),
    is_active: Optional[bool] = None,
    admin_user: User = Depends(require_admin),
    db: Session = Depends(get_db)
):
    """
    List all users (admin only)

    Args:
        limit: Maximum number of users to return
        offset: Number of users to skip
        role_filter: Filter by role (admin, user, viewer)
        is_active: Filter by active status
        admin_user: Current admin user
        db: Database session

    Returns:
        List of users
    """
    query = db.query(User)

    # Apply filters
    if role_filter:
        query = query.filter(User.role == role_filter)

    if is_active is not None:
        query = query.filter(User.is_active == is_active)

    # Order by creation date
    query = query.order_by(User.created_at.desc())

    # Apply pagination
    users = query.offset(offset).limit(limit).all()

    return users


@router.post("/", response_model=UserResponse, status_code=status.HTTP_201_CREATED)
async def create_user(
    user_data: UserCreate,
    admin_user: User = Depends(require_admin),
    db: Session = Depends(get_db)
):
    """
    Create a new user (admin only)

    Args:
        user_data: User creation data
        admin_user: Current admin user
        db: Database session

    Returns:
        Created user
    """
    # Check if username already exists
    existing_user = db.query(User).filter(User.username == user_data.username).first()
    if existing_user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Username '{user_data.username}' already exists"
        )

    # Check if email already exists (if provided)
    if user_data.email:
        existing_email = db.query(User).filter(User.email == user_data.email).first()
        if existing_email:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Email '{user_data.email}' already exists"
            )

    # Create new user
    new_user = User(
        username=user_data.username,
        email=user_data.email,
        password_hash=get_password_hash(user_data.password),
        role=user_data.role,
        max_cores=user_data.max_cores,
        max_jobs=user_data.max_jobs,
        is_active=True
    )

    db.add(new_user)
    db.commit()
    db.refresh(new_user)

    return new_user


@router.get("/me", response_model=UserResponse)
async def get_current_user_info(
    current_user: User = Depends(get_current_active_user)
):
    """
    Get current user information

    Args:
        current_user: Current authenticated user

    Returns:
        Current user details
    """
    return current_user


@router.get("/{user_id}", response_model=UserResponse)
async def get_user(
    user_id: UUID,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Get user by ID

    Args:
        user_id: User UUID
        current_user: Current authenticated user
        db: Database session

    Returns:
        User details
    """
    user = db.query(User).filter(User.id == user_id).first()

    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Users can view their own profile, admins can view anyone
    if user.id != current_user.id and current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to view this user"
        )

    return user


@router.patch("/{user_id}", response_model=UserResponse)
async def update_user(
    user_id: UUID,
    user_update: UserUpdate,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    Update user information

    Args:
        user_id: User UUID
        user_update: User update data
        current_user: Current authenticated user
        db: Database session

    Returns:
        Updated user
    """
    user = db.query(User).filter(User.id == user_id).first()

    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Users can update their own profile (limited fields)
    # Admins can update anyone with all fields
    is_self = user.id == current_user.id
    is_admin = current_user.role == "admin"

    if not (is_self or is_admin):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not authorized to update this user"
        )

    # Regular users can only update email and password
    # Admins can update everything
    update_data = user_update.dict(exclude_unset=True)

    if not is_admin:
        # Non-admins can't change role, is_active, max_cores, max_jobs
        restricted_fields = ["role", "is_active", "max_cores", "max_jobs"]
        for field in restricted_fields:
            if field in update_data:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail=f"Not authorized to update field: {field}"
                )

    # Update password hash if password is provided
    if "password" in update_data:
        password = update_data.pop("password")
        user.password_hash = get_password_hash(password)

    # Update other fields
    for field, value in update_data.items():
        setattr(user, field, value)

    user.updated_at = datetime.utcnow()

    db.commit()
    db.refresh(user)

    return user


@router.delete("/{user_id}")
async def delete_user(
    user_id: UUID,
    admin_user: User = Depends(require_admin),
    db: Session = Depends(get_db)
):
    """
    Delete a user (admin only)

    Note: This performs a soft delete by setting is_active=False

    Args:
        user_id: User UUID
        admin_user: Current admin user
        db: Database session

    Returns:
        Success message
    """
    user = db.query(User).filter(User.id == user_id).first()

    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Prevent self-deletion
    if user.id == admin_user.id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot delete your own account"
        )

    # Soft delete
    user.is_active = False
    user.updated_at = datetime.utcnow()

    db.commit()

    return {"message": f"User {user.username} deactivated successfully"}


@router.post("/{user_id}/activate")
async def activate_user(
    user_id: UUID,
    admin_user: User = Depends(require_admin),
    db: Session = Depends(get_db)
):
    """
    Activate a deactivated user (admin only)

    Args:
        user_id: User UUID
        admin_user: Current admin user
        db: Database session

    Returns:
        Success message
    """
    user = db.query(User).filter(User.id == user_id).first()

    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    if user.is_active:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="User is already active"
        )

    user.is_active = True
    user.updated_at = datetime.utcnow()

    db.commit()

    return {"message": f"User {user.username} activated successfully"}
