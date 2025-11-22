# Multi-stage Dockerfile for FIDASIM core
# Stage 1: Build environment
FROM ubuntu:24.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    gfortran \
    gcc \
    g++ \
    make \
    wget \
    zlib1g-dev \
    libgomp1 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy source code
COPY src/ ./src/
COPY tables/ ./tables/
COPY deps/ ./deps/
COPY makefile ./makefile
COPY VERSION ./VERSION

# Set build environment variables
ENV FC=gfortran
ENV CC=gcc
ENV CXX=g++
ENV USE_OPENMP=y
ENV USE_MPI=n

# Build HDF5 from source (always build in container for compatibility)
RUN echo "Building HDF5 from source for container compatibility" && \
    cd deps && \
    tar -xzf hdf5-1.8.16.tar.gz && \
    cd hdf5-1.8.16 && \
    ./configure --prefix=/build/deps/hdf5 \
                --enable-fortran \
                --enable-shared \
                --with-zlib && \
    make -j$(nproc) && \
    make install

# Set HDF5 paths for FIDASIM build
ENV HDF5_HOME=/build/deps/hdf5
ENV PATH=$HDF5_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$HDF5_HOME/lib:$LD_LIBRARY_PATH

# Build FIDASIM with -O1 for stability in container
RUN echo "Building FIDASIM from source with -O1 optimization" && \
    make clean && \
    sed -i 's/-Ofast/-O1/g' makefile && \
    make src -j$(nproc)

# Build atomic tables (or use pre-built if available)
RUN if [ -f tables/atomic_tables.h5 ]; then \
        echo "Using pre-built atomic tables from tables/atomic_tables.h5"; \
    else \
        echo "Building atomic tables (this may take 30+ minutes)"; \
        make tables; \
    fi

# Stage 2: Runtime environment
FROM ubuntu:24.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    gfortran \
    libgomp1 \
    zlib1g \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for the API (updated for Python 3.12 compatibility)
RUN pip3 install --no-cache-dir --break-system-packages \
    fastapi==0.115.0 \
    uvicorn[standard]==0.32.0 \
    pydantic==2.9.0 \
    h5py==3.11.0 \
    numpy==2.1.0 \
    aiofiles==24.1.0 \
    python-multipart==0.0.12

# Create non-root user with explicit UID to match other containers
# Check if UID 1000 exists and remove it if it does
RUN (getent passwd 1000 && userdel -r $(getent passwd 1000 | cut -d: -f1)) || true && \
    useradd -m -u 1000 -s /bin/bash fidasim

# Set working directory
WORKDIR /app

# Copy FIDASIM executable and required files from builder
COPY --from=builder /build/fidasim ./fidasim
COPY --from=builder /build/tables/atomic_tables.h5 ./tables/atomic_tables.h5
COPY --from=builder /build/deps/hdf5/lib/*.so* /usr/local/lib/

# Update library cache
RUN ldconfig

# Create data directory for shared volume
RUN mkdir -p /data && chown fidasim:fidasim /data

# Copy API application
COPY docker/core/api.py ./api.py
COPY docker/core/entrypoint.sh ./entrypoint.sh

# Make entrypoint executable
RUN chmod +x ./entrypoint.sh ./fidasim

# Set ownership
RUN chown -R fidasim:fidasim /app

# Switch to non-root user
USER fidasim

# Environment variables
# Note: OMP_NUM_THREADS is set by FIDASIM command line argument, not here
ENV FIDASIM_DIR=/app
ENV DATA_DIR=/data
ENV PYTHONUNBUFFERED=1

# Expose FastAPI port
EXPOSE 8001

# Health check
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD python3 -c "import urllib.request; urllib.request.urlopen('http://localhost:8001/health').read()" || exit 1

# Entry point
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["api"]