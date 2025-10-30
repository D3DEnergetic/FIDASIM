/**
 * Job Details Page
 *
 * Shows detailed information about a specific job with real-time updates
 */

import { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import {
  Box,
  Paper,
  Typography,
  Button,
  CircularProgress,
  Alert,
  Chip,
  Grid,
  Card,
  CardContent,
  IconButton,
  Divider,
} from '@mui/material';
import {
  ArrowBack as BackIcon,
  Refresh as RefreshIcon,
  Cancel as CancelIcon,
  Download as DownloadIcon,
  InsertChart as ChartIcon,
} from '@mui/icons-material';
import { useQuery, useQueryClient } from '@tanstack/react-query';
import { api } from '@/services/api';
import { websocketService } from '@/services/websocket';
import type { Job, WebSocketMessage } from '@/types';

export default function JobDetails() {
  const { jobId } = useParams<{ jobId: string }>();
  const navigate = useNavigate();
  const queryClient = useQueryClient();
  const [liveJob, setLiveJob] = useState<Job | null>(null);
  const [logs, setLogs] = useState<string>('');

  // Fetch job details
  const {
    data: job,
    isLoading,
    error,
    refetch,
  } = useQuery<Job>({
    queryKey: ['job', jobId],
    queryFn: () => api.getJob(jobId!),
    enabled: !!jobId,
    refetchInterval: (query) => {
      const data = query.state.data;
      return (data?.status === 'pending' || data?.status === 'running') ? 5000 : false;
    },
  });

  // Fetch job logs (get all logs, not just tail 100)
  const { data: logsData, refetch: refetchLogs } = useQuery<{ logs: string }>({
    queryKey: ['jobLogs', jobId],
    queryFn: () => api.getJobLogs(jobId!, 10000), // Get last 10000 lines
    enabled: !!jobId,
    refetchInterval: () => {
      const jobData = queryClient.getQueryData<Job>(['job', jobId]);
      return (jobData?.status === 'running') ? 5000 : false;
    },
  });

  // Update live job when data changes
  useEffect(() => {
    if (job) {
      setLiveJob(job);
    }
  }, [job]);

  // Update logs
  useEffect(() => {
    if (logsData?.logs) {
      setLogs(logsData.logs);
    }
  }, [logsData]);

  // Subscribe to WebSocket updates
  useEffect(() => {
    const unsubscribe = websocketService.subscribe((message: WebSocketMessage) => {
      if (message.type === 'job_update' && message.job_id === jobId) {
        // Update job status
        setLiveJob((prev) => (prev ? { ...prev, status: message.status! } : null));
        refetch();
      } else if (message.type === 'job_log' && message.job_id === jobId) {
        // Append new log line
        setLogs((prev) => prev + '\n' + message.log);
      }
    });

    return () => unsubscribe();
  }, [jobId, refetch]);

  const handleCancelJob = async () => {
    if (!jobId || !liveJob) return;

    if (confirm(`Are you sure you want to cancel job ${liveJob.run_id}?`)) {
      try {
        await api.cancelJob(jobId);
        refetch();
      } catch (err) {
        console.error('Failed to cancel job:', err);
      }
    }
  };

  const getStatusColor = (status: string): 'success' | 'error' | 'warning' | 'info' | 'default' => {
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

  if (isLoading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', mt: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  if (error || !liveJob) {
    return (
      <Box>
        <Alert severity="error">Job not found or error loading job details.</Alert>
        <Button startIcon={<BackIcon />} onClick={() => navigate('/jobs')} sx={{ mt: 2 }}>
          Back to Jobs
        </Button>
      </Box>
    );
  }

  return (
    <Box>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <IconButton onClick={() => navigate('/jobs')}>
            <BackIcon />
          </IconButton>
          <Typography variant="h4" component="h1">
            {liveJob.run_id}
          </Typography>
          <Chip label={liveJob.status} color={getStatusColor(liveJob.status)} />
        </Box>

        <Box>
          <IconButton onClick={() => refetch()} sx={{ mr: 1 }}>
            <RefreshIcon />
          </IconButton>
          {(liveJob.status === 'completed') && (
            <Button
              variant="contained"
              color="primary"
              startIcon={<ChartIcon />}
              onClick={() => navigate(`/visualization/${jobId}`)}
              sx={{ mr: 1 }}
            >
              Visualize Results
            </Button>
          )}
          {(liveJob.status === 'pending' || liveJob.status === 'running') && (
            <Button
              variant="outlined"
              color="error"
              startIcon={<CancelIcon />}
              onClick={handleCancelJob}
            >
              Cancel Job
            </Button>
          )}
        </Box>
      </Box>

      {/* Job Information Cards */}
      <Grid container spacing={3} sx={{ mb: 3 }}>
        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Job ID
              </Typography>
              <Typography variant="body1" sx={{ fontFamily: 'monospace', wordBreak: 'break-all' }}>
                {liveJob.id}
              </Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Cores Assigned
              </Typography>
              <Typography variant="h5">{liveJob.cores_assigned}</Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Created
              </Typography>
              <Typography variant="body1">
                {liveJob.created_at ? new Date(liveJob.created_at).toLocaleString() : 'N/A'}
              </Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Started
              </Typography>
              <Typography variant="body1">
                {liveJob.started_at ? new Date(liveJob.started_at).toLocaleString() : 'Not started'}
              </Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Completed
              </Typography>
              <Typography variant="body1">
                {liveJob.completed_at ? new Date(liveJob.completed_at).toLocaleString() : liveJob.status === 'running' ? 'Running...' : 'N/A'}
              </Typography>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" variant="body2">
                Duration
              </Typography>
              <Typography variant="body1">
                {liveJob.completed_at && liveJob.started_at
                  ? `${Math.round(
                      (new Date(liveJob.completed_at).getTime() -
                        new Date(liveJob.started_at).getTime()) /
                        60000
                    )} min`
                  : liveJob.started_at
                  ? `${Math.round((Date.now() - new Date(liveJob.started_at).getTime()) / 60000)} min (running)`
                  : 'Not started'}
              </Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      {/* Configuration */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Configuration
        </Typography>
        <Divider sx={{ mb: 2 }} />
        <Grid container spacing={2}>
          {Object.entries(liveJob.config).map(([key, value]) => {
            // Skip null, undefined, or empty values
            if (value === null || value === undefined || value === '') return null;

            // Format value for display
            let displayValue = value;
            if (typeof value === 'object') {
              displayValue = JSON.stringify(value, null, 2);
            } else if (typeof value === 'boolean') {
              displayValue = value ? 'Yes' : 'No';
            } else {
              displayValue = String(value);
            }

            return (
              <Grid item xs={12} md={6} key={key}>
                <Typography variant="body2" color="text.secondary">
                  {key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                </Typography>
                <Typography variant="body1" sx={{ wordBreak: 'break-word' }}>
                  {displayValue}
                </Typography>
              </Grid>
            );
          })}
        </Grid>
      </Paper>

      {/* Error Message */}
      {liveJob.error_message && (
        <Alert severity="error" sx={{ mb: 3 }}>
          <Typography variant="body2" fontWeight="bold">
            Error:
          </Typography>
          <Typography variant="body2">{liveJob.error_message}</Typography>
        </Alert>
      )}

      {/* Output Files */}
      {liveJob.output_files && Array.isArray(liveJob.output_files) && liveJob.output_files.length > 0 && (
        <Paper sx={{ p: 3, mb: 3 }}>
          <Typography variant="h6" gutterBottom>
            Output Files ({liveJob.output_files.length})
          </Typography>
          <Divider sx={{ mb: 2 }} />
          {liveJob.output_files.map((file, index) => (
            <Box key={index} sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 1 }}>
              <Typography variant="body2" fontFamily="monospace">{file}</Typography>
              <IconButton size="small" disabled title="Download functionality coming soon">
                <DownloadIcon />
              </IconButton>
            </Box>
          ))}
        </Paper>
      )}

      {/* Logs */}
      <Paper sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
          <Typography variant="h6">Logs</Typography>
          <Button size="small" onClick={() => refetchLogs()}>
            Refresh Logs
          </Button>
        </Box>
        <Divider sx={{ mb: 2 }} />
        <Box
          sx={{
            bgcolor: '#1e1e1e',
            color: '#d4d4d4',
            p: 2,
            borderRadius: 1,
            fontFamily: 'monospace',
            fontSize: '0.875rem',
            maxHeight: 500,
            overflow: 'auto',
            whiteSpace: 'pre-wrap',
            wordBreak: 'break-word',
          }}
        >
          {logs || 'No logs available yet...'}
        </Box>
      </Paper>
    </Box>
  );
}
