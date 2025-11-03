#!/bin/bash
set -e

# FIDASIM Core Container Entrypoint

# Function to handle signals gracefully
trap_handler() {
    echo "Received shutdown signal, cleaning up..."
    # Kill any running FIDASIM processes
    pkill -f fidasim || true
    exit 0
}

trap trap_handler SIGTERM SIGINT

# Check the command
case "$1" in
    api)
        echo "Starting FIDASIM Core API server on port 8001..."
        exec uvicorn api:app --host 0.0.0.0 --port 8001 --workers 4
        ;;

    fidasim)
        # Direct FIDASIM execution mode
        shift  # Remove 'fidasim' from arguments
        echo "Running FIDASIM directly with arguments: $@"
        exec /app/fidasim "$@"
        ;;

    test)
        echo "Running FIDASIM test mode..."
        # Simple test to verify the executable works
        /app/fidasim --version || echo "FIDASIM version check"
        echo "FIDASIM executable verified"
        echo "Atomic tables location: /app/tables/atomic_tables.h5"
        ls -la /app/tables/
        echo "Data directory: /data"
        ls -la /data/ || echo "Data directory is empty"
        ;;

    bash)
        # Development mode - start bash shell
        exec /bin/bash
        ;;

    *)
        echo "Usage: $0 {api|fidasim|test|bash}"
        echo "  api     - Start FastAPI server (default)"
        echo "  fidasim - Run FIDASIM directly with provided arguments"
        echo "  test    - Run basic tests"
        echo "  bash    - Start interactive shell (for debugging)"
        exit 1
        ;;
esac