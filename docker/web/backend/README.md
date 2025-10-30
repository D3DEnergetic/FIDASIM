# FIDASIM Web Backend

FastAPI-based backend for the FIDASIM web interface.

## Features

- **Authentication**: JWT-based authentication with role-based access control
- **Job Management**: Submit, monitor, and manage FIDASIM simulation jobs
- **User Management**: Admin interface for user management
- **File Handling**: Upload/download input and result files
- **Real-time Updates**: WebSocket support for live job status updates
- **Job Orchestration**: Communicates with fidasim-core and fidasim-preprocessor APIs

## Project Structure

```
app/
├── main.py                      # FastAPI application entry point
├── core/
│   ├── config.py               # Application configuration
│   └── security.py             # Authentication & JWT handling
├── db/
│   └── database.py             # SQLAlchemy database setup
├── models/
│   ├── user.py                 # User database model
│   └── job.py                  # Job database model
├── schemas/
│   ├── auth.py                 # Authentication schemas
│   ├── user.py                 # User validation schemas
│   └── job.py                  # Job validation schemas
├── api/
│   ├── auth.py                 # Authentication endpoints
│   ├── jobs.py                 # Job management endpoints
│   ├── users.py                # User management endpoints
│   ├── files.py                # File upload/download endpoints
│   └── websocket.py            # WebSocket endpoints
└── services/
    ├── job_orchestrator.py     # Job submission & tracking service
    └── websocket_manager.py    # WebSocket connection manager
```

## Setup

1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Configure environment**:
   ```bash
   cp .env.example .env
   # Edit .env with your configuration
   ```

3. **Run the server**:
   ```bash
   # Development
   uvicorn app.main:app --reload --host 0.0.0.0 --port 8000

   # Production
   uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4
   ```

## API Documentation

Once running, access the interactive API documentation:
- Swagger UI: http://localhost:8000/api/docs
- ReDoc: http://localhost:8000/api/redoc

## API Endpoints

### Authentication (`/api/auth`)
- `POST /api/auth/login` - Login and get JWT token
- `GET /api/auth/me` - Get current user info
- `POST /api/auth/logout` - Logout

### Jobs (`/api/jobs`)
- `POST /api/jobs/submit` - Submit new job
- `GET /api/jobs/` - List jobs (with filters)
- `GET /api/jobs/stats` - Get job statistics
- `GET /api/jobs/{id}` - Get job details
- `PATCH /api/jobs/{id}` - Update job
- `DELETE /api/jobs/{id}/cancel` - Cancel job
- `GET /api/jobs/{id}/logs` - Get job logs

### Users (`/api/users`)
- `GET /api/users/` - List all users (admin only)
- `POST /api/users/` - Create user (admin only)
- `GET /api/users/me` - Get current user
- `GET /api/users/{id}` - Get user by ID
- `PATCH /api/users/{id}` - Update user
- `DELETE /api/users/{id}` - Deactivate user (admin only)
- `POST /api/users/{id}/activate` - Activate user (admin only)

### Files (`/api/files`)
- `POST /api/files/upload` - Upload file
- `POST /api/files/upload-multiple` - Upload multiple files
- `GET /api/files/list` - List uploaded files
- `GET /api/files/download/{filename}` - Download file
- `GET /api/files/info/{filename}` - Get file info
- `GET /api/files/storage-usage` - Get storage usage stats
- `DELETE /api/files/{filename}` - Delete file

### WebSocket (`/api/ws`)
- `WS /api/ws/jobs?token=<jwt>` - Real-time job updates
- `WS /api/ws/notifications?token=<jwt>` - Real-time notifications
- `GET /api/ws/stats` - WebSocket connection stats

## Authentication

The API uses JWT (JSON Web Tokens) for authentication:

1. **Login** to get a token:
   ```bash
   curl -X POST http://localhost:8000/api/auth/login \
     -H "Content-Type: application/json" \
     -d '{"username": "admin", "password": "password"}'
   ```

2. **Use the token** in subsequent requests:
   ```bash
   curl http://localhost:8000/api/jobs/ \
     -H "Authorization: Bearer <your_token>"
   ```

## WebSocket Connection

Connect to WebSocket endpoints with JWT token:

```javascript
const token = "your_jwt_token";
const ws = new WebSocket(`ws://localhost:8000/api/ws/jobs?token=${token}`);

ws.onmessage = (event) => {
  const data = JSON.parse(event.data);
  console.log("Received:", data);
};

// Keep connection alive
setInterval(() => ws.send("ping"), 30000);
```

## Development

### Creating a test user

You'll need to create users directly in the database initially. Once you have an admin user, you can use the `/api/users/` endpoint.

### Running tests

```bash
pytest tests/
```

## Environment Variables

See `.env.example` for all available configuration options.

Key variables:
- `SECRET_KEY`: JWT signing key (generate with `openssl rand -hex 32`)
- `DATABASE_URL`: PostgreSQL connection string
- `CORE_API_URL`: URL for fidasim-core service
- `PREPROCESSOR_API_URL`: URL for fidasim-preprocessor service

## Deployment

For production deployment:

1. Set strong `SECRET_KEY`
2. Configure proper CORS origins
3. Use PostgreSQL database (not SQLite)
4. Set up Redis for caching
5. Use reverse proxy (Nginx) with HTTPS
6. Run with multiple workers: `--workers 4`

## Integration with FIDASIM Core

The backend communicates with fidasim-core via HTTP:

```
Web Backend → Job Orchestrator → fidasim-core API
                               → fidasim-preprocessor API
```

The `JobOrchestrator` service handles:
- Job submission to core API
- Status polling
- Log retrieval
- Result file downloads
- Job cancellation
