"""
WebSocket API endpoints

Provides real-time updates for job status, logs, and notifications
"""

from fastapi import APIRouter, WebSocket, WebSocketDisconnect, Depends, Query
from sqlalchemy.orm import Session
from typing import Optional

from app.core.security import decode_access_token
from app.db.database import get_db
from app.models.user import User
from app.services.websocket_manager import websocket_manager

router = APIRouter()


async def get_user_from_token(token: str, db: Session) -> Optional[User]:
    """
    Get user from JWT token for WebSocket authentication

    Args:
        token: JWT token
        db: Database session

    Returns:
        User object or None
    """
    payload = decode_access_token(token)
    if not payload:
        return None

    username = payload.get("sub")
    if not username:
        return None

    user = db.query(User).filter(User.username == username).first()
    if not user or not user.is_active:
        return None

    return user


@router.websocket("/jobs")
async def websocket_job_updates(
    websocket: WebSocket,
    token: str = Query(..., description="JWT authentication token"),
    db: Session = Depends(get_db)
):
    """
    WebSocket endpoint for real-time job updates

    Clients should connect with a valid JWT token as a query parameter:
    ws://host/api/ws/jobs?token=<jwt_token>

    Message types sent to client:
    - job_update: Job status changed
    - job_log: New log line available
    - notification: General notification

    Example messages:
    {
        "type": "job_update",
        "job_id": "uuid",
        "status": "running",
        "timestamp": "2024-01-01T00:00:00",
        "details": {}
    }

    {
        "type": "job_log",
        "job_id": "uuid",
        "log": "Log line content",
        "timestamp": "2024-01-01T00:00:00"
    }

    {
        "type": "notification",
        "notification_type": "info",
        "title": "Notification title",
        "body": "Notification body",
        "timestamp": "2024-01-01T00:00:00"
    }
    """
    # Authenticate user
    user = await get_user_from_token(token, db)
    if not user:
        await websocket.close(code=1008, reason="Invalid authentication token")
        return

    # Connect to WebSocket manager
    await websocket_manager.connect(websocket, str(user.id))

    try:
        # Send welcome message
        await websocket.send_json({
            "type": "connected",
            "message": "Connected to FIDASIM WebSocket",
            "user_id": str(user.id),
            "username": user.username
        })

        # Keep connection alive and handle incoming messages
        while True:
            # Receive messages from client (for future interactivity)
            data = await websocket.receive_text()

            # Echo ping/pong for keepalive
            if data == "ping":
                await websocket.send_json({
                    "type": "pong",
                    "timestamp": "2024-01-01T00:00:00"
                })

    except WebSocketDisconnect:
        # Client disconnected
        await websocket_manager.disconnect(websocket, str(user.id))
    except Exception as e:
        # Handle other errors
        print(f"WebSocket error for user {user.id}: {e}")
        await websocket_manager.disconnect(websocket, str(user.id))
        try:
            await websocket.close(code=1011, reason=f"Server error: {str(e)}")
        except:
            pass


@router.websocket("/notifications")
async def websocket_notifications(
    websocket: WebSocket,
    token: str = Query(..., description="JWT authentication token"),
    db: Session = Depends(get_db)
):
    """
    WebSocket endpoint for general notifications

    Provides system-wide notifications, announcements, and alerts.
    """
    # Authenticate user
    user = await get_user_from_token(token, db)
    if not user:
        await websocket.close(code=1008, reason="Invalid authentication token")
        return

    # Connect to WebSocket manager
    await websocket_manager.connect(websocket, str(user.id))

    try:
        # Send welcome message
        await websocket.send_json({
            "type": "connected",
            "message": "Connected to notifications",
            "user_id": str(user.id)
        })

        # Keep connection alive
        while True:
            data = await websocket.receive_text()

            # Echo ping/pong for keepalive
            if data == "ping":
                await websocket.send_json({"type": "pong"})

    except WebSocketDisconnect:
        await websocket_manager.disconnect(websocket, str(user.id))
    except Exception as e:
        print(f"WebSocket error for user {user.id}: {e}")
        await websocket_manager.disconnect(websocket, str(user.id))


@router.get("/stats")
async def get_websocket_stats():
    """
    Get WebSocket connection statistics

    Returns:
        Connection statistics
    """
    return {
        "total_connections": websocket_manager.get_connection_count(),
        "status": "operational"
    }
