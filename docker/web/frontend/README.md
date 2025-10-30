# FIDASIM Web Frontend

Modern React-based frontend for the FIDASIM simulation platform.

## Features

- **Authentication**: JWT-based login with secure token storage
- **Dashboard**: Real-time job statistics and recent jobs overview
- **Job Management**: Submit, monitor, and manage simulation jobs
- **Live Updates**: WebSocket integration for real-time job status and logs
- **File Management**: Upload and manage input files for simulations
- **User Profile**: Self-service profile management
- **Responsive Design**: Works on desktop, tablet, and mobile devices

## Tech Stack

- **React 18** - UI library
- **TypeScript** - Type safety
- **Vite** - Fast build tool and dev server
- **Material-UI (MUI)** - Component library
- **React Router** - Client-side routing
- **TanStack Query** - Server state management
- **Axios** - HTTP client
- **WebSocket** - Real-time communication

## Project Structure

```
src/
├── main.tsx                  # Application entry point
├── App.tsx                   # Main app component with routing
├── vite-env.d.ts            # TypeScript environment definitions
├── components/
│   ├── Layout.tsx           # Main layout with navigation
│   └── ProtectedRoute.tsx   # Route protection wrapper
├── pages/
│   ├── Login.tsx            # Login page
│   ├── Dashboard.tsx        # Dashboard with job stats
│   ├── Jobs.tsx             # Jobs list page
│   ├── JobDetails.tsx       # Job details with live updates
│   ├── SubmitJob.tsx        # Job submission form
│   ├── Files.tsx            # File upload/management
│   └── Profile.tsx          # User profile page
├── contexts/
│   └── AuthContext.tsx      # Authentication context
├── services/
│   ├── api.ts               # API client
│   └── websocket.ts         # WebSocket service
├── types/
│   └── index.ts             # TypeScript type definitions
└── utils/
    └── (utility functions)
```

## Setup

### Prerequisites

- Node.js 18+ and npm

### Installation

1. **Install dependencies**:
   ```bash
   npm install
   ```

2. **Configure environment**:
   ```bash
   cp .env.example .env
   # Edit .env with your backend API URL
   ```

3. **Run development server**:
   ```bash
   npm run dev
   ```

   The app will be available at http://localhost:3000

### Build for Production

```bash
npm run build
```

The built files will be in the `dist/` directory.

### Preview Production Build

```bash
npm run preview
```

## Environment Variables

See `.env.example` for available configuration options:

- `VITE_API_URL` - Backend API URL (default: http://localhost:8000)
- `VITE_WS_URL` - WebSocket URL (default: ws://localhost:8000)
- `VITE_APP_NAME` - Application name
- `VITE_APP_VERSION` - Application version

## Features Overview

### Authentication

- Secure JWT-based authentication
- Automatic token refresh
- Protected routes
- Role-based access control

### Dashboard

- Real-time job statistics (total, pending, running, completed, failed, cancelled)
- Recent jobs table with live status updates
- Quick access to job submission

### Job Management

- Submit new simulation jobs with customizable parameters
- View all jobs with filtering by status
- Real-time job status updates via WebSocket
- Detailed job view with:
  - Live log streaming
  - Configuration details
  - Error messages
  - Output file downloads
  - Job cancellation

### File Management

- Upload single or multiple files
- View storage usage statistics
- Download files
- Delete files
- Secure user-isolated file storage

### User Profile

- View account information
- Update email
- Change password
- View resource limits (max cores, max jobs)

## Development

### Code Style

- TypeScript strict mode enabled
- ESLint for code linting
- Consistent component structure

### Component Patterns

- Functional components with hooks
- TypeScript interfaces for props
- Material-UI components for consistent design
- React Query for server state
- Context API for global state (auth)

### API Integration

All API calls go through the centralized `api.ts` service:

```typescript
import { api } from '@/services/api';

// Example: Get jobs
const jobs = await api.getJobs();

// Example: Submit job
const job = await api.submitJob(jobData);
```

### WebSocket Integration

Real-time updates are handled by the WebSocket service:

```typescript
import { websocketService } from '@/services/websocket';

// Subscribe to updates
const unsubscribe = websocketService.subscribe((message) => {
  if (message.type === 'job_update') {
    // Handle job update
  }
});

// Cleanup
return () => unsubscribe();
```

## Deployment

### Docker Deployment

A Dockerfile will be provided for containerized deployment.

### Environment Setup

For production:
1. Set `VITE_API_URL` to your production backend URL
2. Set `VITE_WS_URL` to your production WebSocket URL
3. Build with `npm run build`
4. Serve the `dist/` directory with a web server (Nginx, Apache, etc.)

### Nginx Configuration Example

```nginx
server {
    listen 80;
    server_name your-domain.com;

    root /path/to/dist;
    index index.html;

    location / {
        try_files $uri $uri/ /index.html;
    }

    # Proxy API requests to backend
    location /api {
        proxy_pass http://backend:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    # Proxy WebSocket connections
    location /api/ws {
        proxy_pass http://backend:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

## Troubleshooting

### WebSocket Connection Issues

If WebSocket connections fail:
1. Check that `VITE_WS_URL` is correctly configured
2. Ensure backend WebSocket endpoint is accessible
3. Check browser console for connection errors
4. Verify JWT token is being sent correctly

### API Request Failures

If API requests fail:
1. Check that `VITE_API_URL` is correctly configured
2. Verify backend is running and accessible
3. Check browser network tab for detailed error messages
4. Ensure CORS is properly configured on backend

## License

Same as FIDASIM main project.
