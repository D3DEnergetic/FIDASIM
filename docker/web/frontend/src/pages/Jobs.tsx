/**
 * Jobs List Page
 */

import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Box,
  Paper,
  Typography,
  Button,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TablePagination,
  IconButton,
  Chip,
  TextField,
  MenuItem,
  CircularProgress,
  Alert,
} from '@mui/material';
import {
  Visibility as ViewIcon,
  Refresh as RefreshIcon,
  Add as AddIcon,
} from '@mui/icons-material';
import { useQuery } from '@tanstack/react-query';
import { api } from '@/services/api';
import { websocketService } from '@/services/websocket';
import type { Job, JobStatus, WebSocketMessage } from '@/types';

export default function Jobs() {
  const navigate = useNavigate();
  const [page, setPage] = useState(0);
  const [rowsPerPage, setRowsPerPage] = useState(25);
  const [statusFilter, setStatusFilter] = useState<string>('');
  const [liveJobs, setLiveJobs] = useState<Job[]>([]);

  // Fetch jobs
  const {
    data: jobs,
    isLoading,
    error,
    refetch,
  } = useQuery<Job[]>({
    queryKey: ['jobs', statusFilter, rowsPerPage, page * rowsPerPage],
    queryFn: () => api.getJobs(statusFilter || undefined, rowsPerPage, page * rowsPerPage),
    refetchInterval: 10000, // Refetch every 10 seconds
  });

  // Update live jobs when data changes
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
      }
    });

    return () => unsubscribe();
  }, []);

  const handleChangePage = (_event: unknown, newPage: number) => {
    setPage(newPage);
  };

  const handleChangeRowsPerPage = (event: React.ChangeEvent<HTMLInputElement>) => {
    setRowsPerPage(parseInt(event.target.value, 10));
    setPage(0);
  };

  const getStatusColor = (status: JobStatus): 'success' | 'error' | 'warning' | 'info' | 'default' => {
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

  if (error) {
    return <Alert severity="error">Error loading jobs. Please try again.</Alert>;
  }

  return (
    <Box>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
        <Typography variant="h4" component="h1">
          My Jobs
        </Typography>
        <Box>
          <IconButton onClick={() => refetch()} sx={{ mr: 1 }}>
            <RefreshIcon />
          </IconButton>
          <Button variant="contained" startIcon={<AddIcon />} onClick={() => navigate('/submit')}>
            Submit Job
          </Button>
        </Box>
      </Box>

      {/* Filters */}
      <Paper sx={{ p: 2, mb: 3 }}>
        <TextField
          select
          label="Filter by Status"
          value={statusFilter}
          onChange={(e) => {
            setStatusFilter(e.target.value);
            setPage(0);
          }}
          sx={{ minWidth: 200 }}
          size="small"
        >
          <MenuItem value="">All</MenuItem>
          <MenuItem value="pending">Pending</MenuItem>
          <MenuItem value="running">Running</MenuItem>
          <MenuItem value="completed">Completed</MenuItem>
          <MenuItem value="failed">Failed</MenuItem>
          <MenuItem value="cancelled">Cancelled</MenuItem>
        </TextField>
      </Paper>

      {/* Jobs Table */}
      <Paper>
        <TableContainer>
          <Table>
            <TableHead>
              <TableRow>
                <TableCell>Run ID</TableCell>
                <TableCell>Status</TableCell>
                <TableCell>Cores</TableCell>
                <TableCell>Priority</TableCell>
                <TableCell>Created</TableCell>
                <TableCell>Started</TableCell>
                <TableCell>Completed</TableCell>
                <TableCell align="right">Actions</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {liveJobs.length === 0 ? (
                <TableRow>
                  <TableCell colSpan={8} align="center">
                    <Typography color="text.secondary" sx={{ py: 3 }}>
                      {statusFilter ? 'No jobs found with this status' : 'No jobs yet'}
                    </Typography>
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
                      <Chip label={job.status} color={getStatusColor(job.status)} size="small" />
                    </TableCell>
                    <TableCell>{job.cores_assigned}</TableCell>
                    <TableCell>{job.priority}</TableCell>
                    <TableCell>{new Date(job.created_at).toLocaleString()}</TableCell>
                    <TableCell>
                      {job.started_at ? new Date(job.started_at).toLocaleString() : '-'}
                    </TableCell>
                    <TableCell>
                      {job.completed_at ? new Date(job.completed_at).toLocaleString() : '-'}
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

        <TablePagination
          rowsPerPageOptions={[10, 25, 50, 100]}
          component="div"
          count={-1} // Unknown total, using server-side pagination
          rowsPerPage={rowsPerPage}
          page={page}
          onPageChange={handleChangePage}
          onRowsPerPageChange={handleChangeRowsPerPage}
        />
      </Paper>
    </Box>
  );
}
