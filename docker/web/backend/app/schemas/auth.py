"""
Authentication schemas
"""

from __future__ import annotations
from pydantic import BaseModel
from typing import Any


class Token(BaseModel):
    """JWT token response"""
    access_token: str
    token_type: str = "bearer"
    user: Any  # Will be UserResponse, using Any to avoid circular import

    class Config:
        from_attributes = True


class TokenData(BaseModel):
    """Token payload data"""
    username: str


class Login(BaseModel):
    """Login request"""
    username: str
    password: str