/**
 * Dashboard Page
 *
 * Shows job statistics and recent jobs
 */

import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Box,
  Grid,
  Card,
  CardContent,
  Typography,
  Button,
  CircularProgress,
  Alert,
  Chip,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  IconButton,
} from '@mui/material';
import {
  Add as AddIcon,
  Refresh as RefreshIcon,
  CheckCircle as CompletedIcon,
  Error as ErrorIcon,
  HourglassEmpty as PendingIcon,
  PlayArrow as RunningIcon,
  Cancel as CancelledIcon,
  Visibility as ViewIcon,
} from '@mui/icons-material';
import { useQuery } from '@tanstack/react-query';
import { api } from '@/services/api';
import { websocketService } from '@/services/websocket';
import type { Job, JobStats, WebSocketMessage } from '@/types';

export default function Dashboard() {
  const navigate = useNavigate();
  const [liveJobs, setLiveJobs] = useState<Job[]>([]);

  // Fetch job statistics
  const {
    data: stats,
    isLoading: statsLoading,
    error: statsError,
    refetch: refetchStats,
  } = useQuery<JobStats>({
    queryKey: ['jobStats'],
    queryFn: () => api.getJobStats(),
    refetchInterval: 30000, // Refetch every 30 seconds
  });

  // Fetch recent jobs
  const {
    data: jobs,
    isLoading: jobsLoading,
    error: jobsError,
    refetch: refetchJobs,
  } = useQuery<Job[]>({
    queryKey: ['recentJobs'],
    queryFn: () => api.getJobs(undefined, 10, 0),
    refetchInterval: 30000,
  });

  // Update jobs when data changes
  useEffect(() => {
    if (jobs) {
      setLiveJobs(jobs);
    }
  }, [jobs]);

  // Subscribe to WebSocket updates
  useEffect(() => {
    const unsubscribe = websocketService.subscribe((message: WebSocketMessage) => {
      if (message.type === 'job_update' && message.job_id) {
        // Update the job in the list
        setLiveJobs((prevJobs) =>
          prevJobs.map((job) =>
            job.id === message.job_id ? { ...job, status: message.status! } : job
          )
        );

        // Refetch stats to update counters
        refetchStats();
      }
    });

    return () => unsubscribe();
  }, [refetchStats]);

  const handleRefresh = () => {
    refetchStats();
    refetchJobs();
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CompletedIcon color="success" />;
      case 'failed':
        return <ErrorIcon color="error" />;
      case 'running':
        return <RunningIcon color="info" />;
      case 'pending':
        return <PendingIcon color="warning" />;
      case 'cancelled':
        return <CancelledIcon color="disabled" />;
      default:
        return null;
    }
  };

  const getStatusColor = (
    status: string
  ): 'default' | 'success' | 'error' | 'warning' | 'info' => {
    switch (status) {
      case 'completed':
        return 'success';
      case 'failed':
        return 'error';
      case 'running':
        return 'info';
      case 'pending':
        return 'warning';
      default:
        return 'default';
    }
  };

  if (statsLoading || jobsLoading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', mt: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  if (statsError || jobsError) {
    return (
      <Alert severity="error">
        Error loading dashboard data. Please try refreshing the page.
      </Alert>
    );
  }

  return (
    <Box>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
        <Typography variant="h4" component="h1">
          Dashboard
        </Typography>
        <Box>
          <IconButton onClick={handleRefresh} sx={{ mr: 1 }}>
            <RefreshIcon />
          </IconButton>
          <Button
            variant="contained"
            startIcon={<AddIcon />}
            onClick={() => navigate('/submit')}
          >
            Submit Job
          </Button>
        </Box>
      </Box>

      {/* Job Statistics Cards */}
      <Grid container spacing={3} sx={{ mb: 4 }}>
        <Grid item xs={12} sm={6} md={2}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" gutterBottom>
                Total Jobs
              </Typography>
              <Typography variant="h4">{stats?.total || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={2}>
          <Card sx={{ borderLeft: 4, borderColor: 'warning.main' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <PendingIcon color="warning" />
                <Typography color="text.secondary">Pending</Typography>
              </Box>
              <Typography variant="h4">{stats?.pending || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={2}>
          <Card sx={{ borderLeft: 4, borderColor: 'info.main' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <RunningIcon color="info" />
                <Typography color="text.secondary">Running</Typography>
              </Box>
              <Typography variant="h4">{stats?.running || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={2}>
          <Card sx={{ borderLeft: 4, borderColor: 'success.main' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <CompletedIcon color="success" />
                <Typography color="text.secondary">Completed</Typography>
              </Box>
              <Typography variant="h4">{stats?.completed || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={2}>
          <Card sx={{ borderLeft: 4, borderColor: 'error.main' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <ErrorIcon color="error" />
                <Typography color="text.secondary">Failed</Typography>
              </Box>
              <Typography variant="h4">{stats?.failed || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={2}>
          <Card>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <CancelledIcon />
                <Typography color="text.secondary">Cancelled</Typography>
              </Box>
              <Typography variant="h4">{stats?.cancelled || 0}</Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      {/* Recent Jobs Table */}
      <Paper>
        <Box sx={{ p: 2, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <Typography variant="h6">Recent Jobs</Typography>
          <Button size="small" onClick={() => navigate('/jobs')}>
            View All
          </Button>
        </Box>

        <TableContainer>
          <Table>
            <TableHead>
              <TableRow>
                <TableCell>Run ID</TableCell>
                <TableCell>Status</TableCell>
                <TableCell>Cores</TableCell>
                <TableCell>Created</TableCell>
                <TableCell>Duration</TableCell>
                <TableCell align="right">Actions</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {liveJobs.length === 0 ? (
                <TableRow>
                  <TableCell colSpan={6} align="center">
                    <Typography color="text.secondary">No jobs yet</Typography>
                  </TableCell>
                </TableRow>
              ) : (
                liveJobs.map((job) => (
                  <TableRow key={job.id} hover>
                    <TableCell>
                      <Typography variant="body2" fontWeight="medium">
                        {job.run_id}
                      </Typography>
                    </TableCell>
                    <TableCell>
                      <Chip
                        icon={getStatusIcon(job.status) || undefined}
                        label={job.status}
                        color={getStatusColor(job.status)}
                        size="small"
                      />
                    </TableCell>
                    <TableCell>{job.cores_assigned}</TableCell>
                    <TableCell>
                      {new Date(job.created_at).toLocaleDateString()}
                    </TableCell>
                    <TableCell>
                      {job.completed_at && job.started_at
                        ? `${Math.round(
                            (new Date(job.completed_at).getTime() -
                              new Date(job.started_at).getTime()) /
                              60000
                          )} min`
                        : job.started_at
                        ? 'Running...'
                        : '-'}
                    </TableCell>
                    <TableCell align="right">
                      <IconButton size="small" onClick={() => navigate(`/jobs/${job.id}`)}>
                        <ViewIcon />
                      </IconButton>
                    </TableCell>
                  </TableRow>
                ))
              )}
            </TableBody>
          </Table>
        </TableContainer>
      </Paper>
    </Box>
  );
}
