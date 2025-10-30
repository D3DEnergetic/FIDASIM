"""
Job model
"""

from sqlalchemy import Column, String, Integer, Float, BigInteger, Text, DateTime, JSON, ForeignKey
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.sql import func
from sqlalchemy.orm import relationship
import uuid

from app.db.database import Base


class Job(Base):
    """Job model matching existing PostgreSQL schema"""

    __tablename__ = "jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id", ondelete="CASCADE"), nullable=True)
    run_id = Column(String(255), unique=True, nullable=False, index=True)
    status = Column(String(20), default="pending", index=True)  # pending, running, completed, failed, cancelled
    cores_assigned = Column(Integer, default=4)
    priority = Column(Integer, default=5)

    # Input configuration
    input_files = Column(JSON, default=dict)
    config = Column(JSON, nullable=False)
    namelist_content = Column(Text, nullable=True)

    # Execution details
    container_id = Column(String(100), nullable=True)
    node_assigned = Column(String(100), nullable=True)
    pid = Column(Integer, nullable=True)

    # Timing
    created_at = Column(DateTime(timezone=True), server_default=func.now(), index=True)
    queued_at = Column(DateTime(timezone=True), nullable=True)
    started_at = Column(DateTime(timezone=True), nullable=True)
    completed_at = Column(DateTime(timezone=True), nullable=True)

    # Results
    output_files = Column(JSON, default=list)
    error_message = Column(Text, nullable=True)
    exit_code = Column(Integer, nullable=True)

    # Performance metrics
    cpu_time_seconds = Column(Float, nullable=True)
    memory_peak_mb = Column(Float, nullable=True)
    particles_simulated = Column(BigInteger, nullable=True)

    # Metadata (using metadata_ to avoid SQLAlchemy reserved word)
    metadata_ = Column("metadata", JSON, default=dict)

    # Relationship to user
    user = relationship("User", backref="jobs")