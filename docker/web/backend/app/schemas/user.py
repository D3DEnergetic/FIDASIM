"""
User schemas for API validation
"""

from pydantic import BaseModel, EmailStr, Field
from typing import Optional
from datetime import datetime
from uuid import UUID


class UserBase(BaseModel):
    """Base user schema"""
    username: str = Field(..., min_length=3, max_length=100)
    email: Optional[EmailStr] = None
    role: str = Field(default="user", pattern="^(admin|user|viewer)$")
    max_cores: int = Field(default=32, ge=1, le=128)
    max_jobs: int = Field(default=10, ge=1, le=1000)


class UserCreate(UserBase):
    """Schema for creating a user"""
    password: str = Field(..., min_length=8)


class UserUpdate(BaseModel):
    """Schema for updating a user"""
    email: Optional[EmailStr] = None
    password: Optional[str] = Field(None, min_length=8)
    role: Optional[str] = Field(None, pattern="^(admin|user|viewer)$")
    is_active: Optional[bool] = None
    max_cores: Optional[int] = Field(None, ge=1, le=128)
    max_jobs: Optional[int] = Field(None, ge=1, le=1000)


class UserInDB(UserBase):
    """User schema with all database fields"""
    id: UUID
    is_active: bool
    created_at: datetime
    updated_at: datetime
    last_login: Optional[datetime] = None

    class Config:
        from_attributes = True


class UserResponse(UserInDB):
    """User response schema (excludes sensitive data)"""
    pass