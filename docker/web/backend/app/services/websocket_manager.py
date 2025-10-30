"""
WebSocket Manager Service

Manages WebSocket connections and broadcasts real-time updates to clients.
"""

import json
import asyncio
from typing import Dict, Set, Any, Optional
from uuid import UUID
from fastapi import WebSocket
from datetime import datetime


class ConnectionManager:
    """
    Manages WebSocket connections and message broadcasting
    """

    def __init__(self):
        # Dictionary mapping user_id to set of WebSocket connections
        self.active_connections: Dict[str, Set[WebSocket]] = {}
        # Lock for thread-safe operations
        self._lock = asyncio.Lock()

    async def connect(self, websocket: WebSocket, user_id: str):
        """
        Accept a new WebSocket connection

        Args:
            websocket: WebSocket connection
            user_id: User ID string
        """
        await websocket.accept()

        async with self._lock:
            if user_id not in self.active_connections:
                self.active_connections[user_id] = set()
            self.active_connections[user_id].add(websocket)

        print(f"WebSocket connected for user {user_id}. Total connections: {self.get_connection_count()}")

    async def disconnect(self, websocket: WebSocket, user_id: str):
        """
        Remove a WebSocket connection

        Args:
            websocket: WebSocket connection
            user_id: User ID string
        """
        async with self._lock:
            if user_id in self.active_connections:
                self.active_connections[user_id].discard(websocket)
                # Clean up empty sets
                if not self.active_connections[user_id]:
                    del self.active_connections[user_id]

        print(f"WebSocket disconnected for user {user_id}. Total connections: {self.get_connection_count()}")

    async def send_personal_message(self, message: Dict[str, Any], user_id: str):
        """
        Send a message to all connections of a specific user

        Args:
            message: Message dictionary
            user_id: User ID string
        """
        if user_id not in self.active_connections:
            return

        # Create list of connections to avoid modifying set during iteration
        connections = list(self.active_connections[user_id])
        disconnected = []

        for connection in connections:
            try:
                await connection.send_json(message)
            except Exception as e:
                print(f"Error sending message to user {user_id}: {e}")
                disconnected.append(connection)

        # Clean up disconnected connections
        if disconnected:
            async with self._lock:
                for connection in disconnected:
                    if user_id in self.active_connections:
                        self.active_connections[user_id].discard(connection)

    async def broadcast_to_all(self, message: Dict[str, Any]):
        """
        Broadcast a message to all connected clients

        Args:
            message: Message dictionary
        """
        for user_id in list(self.active_connections.keys()):
            await self.send_personal_message(message, user_id)

    async def broadcast_job_update(self, job_id: UUID, user_id: str, status: str, details: Optional[Dict[str, Any]] = None):
        """
        Broadcast a job status update to the job owner

        Args:
            job_id: Job UUID
            user_id: User ID who owns the job
            status: New job status
            details: Additional job details
        """
        message = {
            "type": "job_update",
            "job_id": str(job_id),
            "status": status,
            "timestamp": datetime.utcnow().isoformat(),
            "details": details or {}
        }

        await self.send_personal_message(message, str(user_id))

    async def broadcast_job_log(self, job_id: UUID, user_id: str, log_line: str):
        """
        Broadcast a job log line to the job owner

        Args:
            job_id: Job UUID
            user_id: User ID who owns the job
            log_line: Log line to send
        """
        message = {
            "type": "job_log",
            "job_id": str(job_id),
            "log": log_line,
            "timestamp": datetime.utcnow().isoformat()
        }

        await self.send_personal_message(message, str(user_id))

    async def send_notification(self, user_id: str, notification_type: str, title: str, body: str):
        """
        Send a notification to a user

        Args:
            user_id: User ID string
            notification_type: Type of notification (info, warning, error, success)
            title: Notification title
            body: Notification body
        """
        message = {
            "type": "notification",
            "notification_type": notification_type,
            "title": title,
            "body": body,
            "timestamp": datetime.utcnow().isoformat()
        }

        await self.send_personal_message(message, user_id)

    def get_connection_count(self) -> int:
        """
        Get total number of active connections

        Returns:
            Number of active WebSocket connections
        """
        return sum(len(connections) for connections in self.active_connections.values())

    def get_user_connection_count(self, user_id: str) -> int:
        """
        Get number of connections for a specific user

        Args:
            user_id: User ID string

        Returns:
            Number of connections for the user
        """
        return len(self.active_connections.get(user_id, set()))

    def is_user_connected(self, user_id: str) -> bool:
        """
        Check if a user has any active connections

        Args:
            user_id: User ID string

        Returns:
            True if user has at least one connection
        """
        return user_id in self.active_connections and len(self.active_connections[user_id]) > 0


# Global WebSocket manager instance
websocket_manager = ConnectionManager()
