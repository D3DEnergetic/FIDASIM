/**
 * Job Submission Page
 */

import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Box,
  Paper,
  Typography,
  TextField,
  Button,
  Alert,
  CircularProgress,
  Grid,
  FormControlLabel,
  Switch,
  Slider,
} from '@mui/material';
import { Send as SendIcon } from '@mui/icons-material';
import { api } from '@/services/api';
import type { JobSubmit } from '@/types';

export default function SubmitJob() {
  const navigate = useNavigate();
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState(false);

  const [formData, setFormData] = useState({
    run_id: '',
    input_file: '',
    device: 'fidasim_device',
    result_dir: '',
    tables_file: '/app/tables/atomic_tables.h5',
    passive_only: false,
    cores: 4,
  });

  const handleChange = (field: string, value: any) => {
    setFormData((prev) => ({ ...prev, [field]: value }));
    setError('');
    setSuccess(false);
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setSuccess(false);
    setIsSubmitting(true);

    try {
      const jobData: JobSubmit = {
        run_id: formData.run_id,
        config: {
          input_file: formData.input_file,
          device: formData.device,
          result_dir: formData.result_dir || `/data/results/${formData.run_id}`,
          tables_file: formData.tables_file,
          passive_only: formData.passive_only,
          cores: formData.cores,
        },
        priority: 0,
      };

      const job = await api.submitJob(jobData);
      setSuccess(true);

      // Redirect to job details after 2 seconds
      setTimeout(() => {
        navigate(`/jobs/${job.id}`);
      }, 2000);
    } catch (err: any) {
      console.error('Job submission failed:', err);
      setError(err.response?.data?.detail || 'Failed to submit job. Please try again.');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <Box>
      <Typography variant="h4" component="h1" gutterBottom>
        Submit New Job
      </Typography>

      <Paper sx={{ p: 3, mt: 3 }}>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {success && (
          <Alert severity="success" sx={{ mb: 2 }}>
            Job submitted successfully! Redirecting to job details...
          </Alert>
        )}

        <Box component="form" onSubmit={handleSubmit}>
          <Grid container spacing={3}>
            {/* Run ID */}
            <Grid item xs={12}>
              <TextField
                required
                fullWidth
                label="Run ID"
                helperText="Unique identifier for this simulation run"
                value={formData.run_id}
                onChange={(e) => handleChange('run_id', e.target.value)}
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Input File */}
            <Grid item xs={12}>
              <TextField
                required
                fullWidth
                label="Input File"
                helperText="Path to the input file (e.g., inputs.dat)"
                value={formData.input_file}
                onChange={(e) => handleChange('input_file', e.target.value)}
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Device */}
            <Grid item xs={12} md={6}>
              <TextField
                fullWidth
                label="Device Name"
                value={formData.device}
                onChange={(e) => handleChange('device', e.target.value)}
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Result Directory */}
            <Grid item xs={12} md={6}>
              <TextField
                fullWidth
                label="Result Directory"
                helperText="Leave empty to use default"
                value={formData.result_dir}
                onChange={(e) => handleChange('result_dir', e.target.value)}
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Atomic Tables File */}
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Atomic Tables File"
                value={formData.tables_file}
                onChange={(e) => handleChange('tables_file', e.target.value)}
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Cores Slider */}
            <Grid item xs={12}>
              <Typography gutterBottom>Number of Cores: {formData.cores}</Typography>
              <Slider
                value={formData.cores}
                onChange={(_, value) => handleChange('cores', value)}
                min={1}
                max={32}
                step={1}
                marks={[
                  { value: 1, label: '1' },
                  { value: 8, label: '8' },
                  { value: 16, label: '16' },
                  { value: 32, label: '32' },
                ]}
                valueLabelDisplay="auto"
                disabled={isSubmitting || success}
              />
            </Grid>

            {/* Passive Only */}
            <Grid item xs={12}>
              <FormControlLabel
                control={
                  <Switch
                    checked={formData.passive_only}
                    onChange={(e) => handleChange('passive_only', e.target.checked)}
                    disabled={isSubmitting || success}
                  />
                }
                label="Passive diagnostics only (no beam)"
              />
            </Grid>

            {/* Submit Button */}
            <Grid item xs={12}>
              <Button
                type="submit"
                variant="contained"
                size="large"
                fullWidth
                disabled={isSubmitting || success}
                startIcon={isSubmitting ? <CircularProgress size={20} /> : <SendIcon />}
              >
                {isSubmitting ? 'Submitting...' : 'Submit Job'}
              </Button>
            </Grid>
          </Grid>
        </Box>
      </Paper>

      <Box sx={{ mt: 2 }}>
        <Typography variant="body2" color="text.secondary">
          <strong>Note:</strong> Make sure you have uploaded your input file before submitting the
          job. You can upload files from the "Upload Files" page.
        </Typography>
      </Box>
    </Box>
  );
}
