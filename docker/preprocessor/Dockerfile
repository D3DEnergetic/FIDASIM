# Dockerfile for FIDASIM Preprocessor Service
FROM python:3.10-slim

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    gfortran \
    make \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user with explicit UID to match other containers
# Check if UID 1000 exists and remove it if it does
RUN (getent passwd 1000 && userdel -r $(getent passwd 1000 | cut -d: -f1)) || true && \
    useradd -m -u 1000 -s /bin/bash preprocessor

# Set working directory
WORKDIR /app

# Copy Python requirements first (for better caching)
COPY docker/preprocessor/requirements.txt ./requirements.txt

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy FIDASIM Python libraries
COPY lib/python/ ./lib/python/

# Copy EFIT dependencies
COPY deps/efit/ ./deps/efit/

# Copy the API application
COPY docker/preprocessor/api.py ./api.py

# Create data directory for shared volume
RUN mkdir -p /data && chown preprocessor:preprocessor /data

# Set Python path to include FIDASIM libraries
ENV PYTHONPATH=/app/lib/python
ENV PYTHONUNBUFFERED=1
ENV DATA_DIR=/data

# Set ownership
RUN chown -R preprocessor:preprocessor /app

# Switch to non-root user
USER preprocessor

# Expose FastAPI port
EXPOSE 8002

# Health check
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8002/health').read()" || exit 1

# Run the API
CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8002"]