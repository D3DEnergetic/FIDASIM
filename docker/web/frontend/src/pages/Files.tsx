/**
 * Files Upload/Management Page
 */

import React, { useState } from 'react';
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
  IconButton,
  Alert,
  CircularProgress,
  Card,
  CardContent,
  Grid,
  LinearProgress,
} from '@mui/material';
import {
  CloudUpload as UploadIcon,
  Delete as DeleteIcon,
  Download as DownloadIcon,
  Refresh as RefreshIcon,
  Info as InfoIcon,
} from '@mui/icons-material';
import { useQuery } from '@tanstack/react-query';
import { api } from '@/services/api';
import type { FileInfo, StorageUsage } from '@/types';

export default function Files() {
  const [uploadingFile, setUploadingFile] = useState(false);
  const [uploadError, setUploadError] = useState('');
  const [uploadSuccess, setUploadSuccess] = useState('');

  // Fetch file list
  const {
    data: filesData,
    isLoading: filesLoading,
    error: filesError,
    refetch: refetchFiles,
  } = useQuery<{ files: FileInfo[] }>({
    queryKey: ['files'],
    queryFn: () => api.listFiles(),
  });

  // Fetch storage usage
  const {
    data: storageData,
    isLoading: storageLoading,
    refetch: refetchStorage,
  } = useQuery<StorageUsage>({
    queryKey: ['storageUsage'],
    queryFn: () => api.getStorageUsage(),
  });

  const handleFileUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = event.target.files;
    if (!files || files.length === 0) return;

    setUploadError('');
    setUploadSuccess('');
    setUploadingFile(true);

    try {
      if (files.length === 1) {
        await api.uploadFile(files[0]);
        setUploadSuccess(`File "${files[0].name}" uploaded successfully`);
      } else {
        await api.uploadMultipleFiles(Array.from(files));
        setUploadSuccess(`${files.length} files uploaded successfully`);
      }

      refetchFiles();
      refetchStorage();
    } catch (err: any) {
      console.error('File upload failed:', err);
      setUploadError(err.response?.data?.detail || 'Failed to upload file(s)');
    } finally {
      setUploadingFile(false);
      event.target.value = ''; // Reset input
    }
  };

  const handleDeleteFile = async (filename: string) => {
    if (!confirm(`Are you sure you want to delete "${filename}"?`)) return;

    try {
      await api.deleteFile(filename);
      setUploadSuccess(`File "${filename}" deleted successfully`);
      refetchFiles();
      refetchStorage();
    } catch (err: any) {
      console.error('File deletion failed:', err);
      setUploadError(err.response?.data?.detail || 'Failed to delete file');
    }
  };

  const handleDownloadFile = async (filename: string) => {
    try {
      const blob = await api.downloadFile(filename);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      document.body.removeChild(a);
    } catch (err: any) {
      console.error('File download failed:', err);
      setUploadError('Failed to download file');
    }
  };

  const storagePercentage = storageData
    ? (storageData.total_size / storageData.max_size) * 100
    : 0;

  return (
    <Box>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
        <Typography variant="h4" component="h1">
          File Management
        </Typography>
        <Box>
          <IconButton onClick={() => { refetchFiles(); refetchStorage(); }} sx={{ mr: 1 }}>
            <RefreshIcon />
          </IconButton>
          <Button
            variant="contained"
            component="label"
            startIcon={uploadingFile ? <CircularProgress size={20} /> : <UploadIcon />}
            disabled={uploadingFile}
          >
            {uploadingFile ? 'Uploading...' : 'Upload Files'}
            <input
              type="file"
              hidden
              multiple
              onChange={handleFileUpload}
              disabled={uploadingFile}
            />
          </Button>
        </Box>
      </Box>

      {uploadError && (
        <Alert severity="error" onClose={() => setUploadError('')} sx={{ mb: 2 }}>
          {uploadError}
        </Alert>
      )}

      {uploadSuccess && (
        <Alert severity="success" onClose={() => setUploadSuccess('')} sx={{ mb: 2 }}>
          {uploadSuccess}
        </Alert>
      )}

      {/* Storage Usage */}
      {!storageLoading && storageData && (
        <Grid container spacing={3} sx={{ mb: 3 }}>
          <Grid item xs={12} md={4}>
            <Card>
              <CardContent>
                <Typography color="text.secondary" gutterBottom>
                  Storage Used
                </Typography>
                <Typography variant="h5">{storageData.total_size_human}</Typography>
                <Typography variant="caption" color="text.secondary">
                  of {storageData.max_size_human}
                </Typography>
                <LinearProgress
                  variant="determinate"
                  value={Math.min(storagePercentage, 100)}
                  sx={{ mt: 1 }}
                  color={storagePercentage > 80 ? 'error' : storagePercentage > 50 ? 'warning' : 'primary'}
                />
              </CardContent>
            </Card>
          </Grid>

          <Grid item xs={12} md={4}>
            <Card>
              <CardContent>
                <Typography color="text.secondary" gutterBottom>
                  Total Files
                </Typography>
                <Typography variant="h5">{storageData.file_count}</Typography>
              </CardContent>
            </Card>
          </Grid>
        </Grid>
      )}

      {/* Files Table */}
      <Paper>
        {filesLoading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', p: 4 }}>
            <CircularProgress />
          </Box>
        ) : filesError ? (
          <Alert severity="error" sx={{ m: 2 }}>
            Error loading files
          </Alert>
        ) : (
          <TableContainer>
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>Filename</TableCell>
                  <TableCell>Size</TableCell>
                  <TableCell>Modified</TableCell>
                  <TableCell align="right">Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {!filesData || filesData.files.length === 0 ? (
                  <TableRow>
                    <TableCell colSpan={4} align="center">
                      <Typography color="text.secondary" sx={{ py: 3 }}>
                        No files uploaded yet
                      </Typography>
                    </TableCell>
                  </TableRow>
                ) : (
                  filesData.files.map((file) => (
                    <TableRow key={file.filename} hover>
                      <TableCell>
                        <Typography variant="body2" fontWeight="medium">
                          {file.filename}
                        </Typography>
                      </TableCell>
                      <TableCell>
                        {file.size > 1024 * 1024
                          ? `${(file.size / (1024 * 1024)).toFixed(2)} MB`
                          : `${(file.size / 1024).toFixed(2)} KB`}
                      </TableCell>
                      <TableCell>
                        {file.modified ? new Date(file.modified).toLocaleString() : '-'}
                      </TableCell>
                      <TableCell align="right">
                        <IconButton
                          size="small"
                          onClick={() => handleDownloadFile(file.filename)}
                          title="Download"
                        >
                          <DownloadIcon />
                        </IconButton>
                        <IconButton
                          size="small"
                          color="error"
                          onClick={() => handleDeleteFile(file.filename)}
                          title="Delete"
                        >
                          <DeleteIcon />
                        </IconButton>
                      </TableCell>
                    </TableRow>
                  ))
                )}
              </TableBody>
            </Table>
          </TableContainer>
        )}
      </Paper>

      <Box sx={{ mt: 2 }}>
        <Alert severity="info" icon={<InfoIcon />}>
          <Typography variant="body2">
            <strong>Tip:</strong> Upload your FIDASIM input files here before submitting jobs. Supported file types include .dat, .h5, and .txt files.
          </Typography>
        </Alert>
      </Box>
    </Box>
  );
}
