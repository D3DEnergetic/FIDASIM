"""
Authentication API endpoints
"""

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import HTTPBasic
from sqlalchemy.orm import Session
from datetime import timedelta

from app.core.config import settings
from app.core.security import (
    authenticate_user,
    create_access_token,
    get_current_user_basic,
    get_current_active_user
)
from app.db.database import get_db
from app.schemas.auth import Token, Login
from app.schemas.user import UserResponse
from app.models.user import User
from sqlalchemy.sql import func

router = APIRouter()
security = HTTPBasic()


@router.post("/login", response_model=Token)
async def login(login_data: Login, db: Session = Depends(get_db)):
    """
    Login with username and password to get JWT token
    """
    user = authenticate_user(db, login_data.username, login_data.password)

    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Update last login
    user.last_login = func.now()
    db.commit()

    # Create access token
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": user.username}, expires_delta=access_token_expires
    )

    # Convert user to response model
    user_response = UserResponse.model_validate(user)

    return {"access_token": access_token, "token_type": "bearer", "user": user_response}


@router.get("/me", response_model=UserResponse)
async def get_current_user_info(current_user: User = Depends(get_current_active_user)):
    """
    Get current user information
    """
    return current_user


@router.post("/logout")
async def logout(current_user: User = Depends(get_current_active_user)):
    """
    Logout (client should discard token)
    """
    return {"message": "Successfully logged out"}