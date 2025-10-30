"""
Configuration management for FIDASIM Web API
"""

from pydantic_settings import BaseSettings
from typing import Optional


class Settings(BaseSettings):
    """Application settings"""

    # Application
    APP_NAME: str = "FIDASIM Web"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False

    # Security
    SECRET_KEY: str
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24  # 24 hours

    # Database
    DATABASE_URL: str = "postgresql://fidasim:fidasim_secure_password_change_me@postgres:5432/fidasim"

    # Redis
    REDIS_URL: str = "redis://redis:6379"

    # External Services
    CORE_API_URL: str = "http://fidasim-core:8001"
    PREPROCESSOR_API_URL: str = "http://fidasim-preprocessor:8002"

    # CORS
    CORS_ORIGINS: list[str] = ["http://localhost:3000", "http://localhost:8000"]

    # File uploads
    MAX_UPLOAD_SIZE: int = 1024 * 1024 * 1000  # 1GB
    UPLOAD_DIR: str = "/data/uploads"

    # Job configuration
    MAX_CORES: int = 32
    DEFAULT_CORES: int = 4
    MAX_CONCURRENT_JOBS: int = 10

    class Config:
        env_file = ".env"
        case_sensitive = True


settings = Settings()