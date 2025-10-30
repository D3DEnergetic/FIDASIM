/**
 * Visualization page for exploring FIDASIM output data
 */

import React, { useState, useEffect } from 'react';
import {
  Container,
  Grid,
  Typography,
  Paper,
  Box,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  SelectChangeEvent,
  Alert,
  Breadcrumbs,
  Link,
  Chip,
  CircularProgress,
  Button,
  IconButton,
  Tooltip,
} from '@mui/material';
import { useQuery, useQueryClient } from '@tanstack/react-query';
import { useNavigate, useParams } from 'react-router-dom';
import visualizationService, {
  DatasetInfo,
  PlotData,
  VisualizationFile,
} from '../services/visualization';
import InteractivePlot from '../components/InteractivePlot';
import CoordinateSelector from '../components/CoordinateSelector';
import ErrorBoundary from '../components/ErrorBoundary';
import RefreshIcon from '@mui/icons-material/Refresh';
import InsertChartIcon from '@mui/icons-material/InsertChart';
import FolderOpenIcon from '@mui/icons-material/FolderOpen';
import DownloadIcon from '@mui/icons-material/Download';

const Visualization: React.FC = () => {
  const navigate = useNavigate();
  const { jobId } = useParams<{ jobId?: string }>();
  const queryClient = useQueryClient();

  // State
  const [selectedFile, setSelectedFile] = useState<string>('');
  const [selectedVariable, setSelectedVariable] = useState<string>('');
  const [xDimension, setXDimension] = useState<string>('');
  const [yDimension, setYDimension] = useState<string | null>(null);
  const [fixedCoords, setFixedCoords] = useState<Record<string, number>>({});
  const [integrateDims, setIntegrateDims] = useState<string[]>([]);
  const [plotData, setPlotData] = useState<PlotData | null>(null);
  const [plotLoading, setPlotLoading] = useState<boolean>(false);
  const [plotError, setPlotError] = useState<string>('');

  // Fetch available files for the job
  const {
    data: files,
    isLoading: filesLoading,
    error: filesError,
  } = useQuery<VisualizationFile[]>({
    queryKey: ['visualization-files', jobId],
    queryFn: () => visualizationService.listFiles(jobId || ''),
    enabled: !!jobId,
  });

  // Fetch dataset info when a file is selected
  const {
    data: datasetInfo,
    isLoading: datasetLoading,
    error: datasetError,
    refetch: refetchDatasetInfo,
  } = useQuery<DatasetInfo>({
    queryKey: ['dataset-info', selectedFile],
    queryFn: () => visualizationService.getDatasetInfo(selectedFile),
    enabled: !!selectedFile,
  });

  // Reset state when file changes
  useEffect(() => {
    // Clear all dependent state when file changes
    setSelectedVariable('');
    setXDimension('');
    setYDimension(null);
    setFixedCoords({});
    setIntegrateDims([]);
    setPlotData(null);
    setPlotError('');
  }, [selectedFile]);

  // Handle file selection
  const handleFileSelect = (event: SelectChangeEvent<string>) => {
    const filePath = event.target.value;
    if (filePath !== selectedFile) {
      setSelectedFile(filePath);
    }
  };

  // Handle visualization
  const handleVisualize = async () => {
    if (!selectedFile || !selectedVariable || !xDimension) {
      return;
    }

    setPlotLoading(true);
    setPlotError('');

    try {
      const data = await visualizationService.getPlotData(
        selectedFile,
        selectedVariable,
        xDimension,
        yDimension || undefined,
        fixedCoords,
        integrateDims.length > 0 ? integrateDims : undefined
      );
      setPlotData(data);
    } catch (error: any) {
      setPlotError(error.response?.data?.detail || error.message || 'Failed to load plot data');
    } finally {
      setPlotLoading(false);
    }
  };

  // Handle dimension changes from CoordinateSelector
  const handleDimensionsChange = (xDim: string, yDim: string | null) => {
    setXDimension(xDim);
    setYDimension(yDim);
  };

  // Handle fixed coordinates change
  const handleFixedCoordsChange = (coords: Record<string, number>) => {
    setFixedCoords(coords);
  };

  // Handle integrate dimensions change
  const handleIntegrateDimsChange = (dims: string[]) => {
    setIntegrateDims(dims);
  };

  // Handle cache clear
  const handleClearCache = async () => {
    await visualizationService.clearCache(selectedFile || undefined);
    queryClient.invalidateQueries({ queryKey: ['dataset-info'] });
    refetchDatasetInfo();
  };

  // Download current plot as image
  const handleDownloadPlot = () => {
    const plotElement = document.querySelector('.js-plotly-plot') as any;
    if (plotElement && window.Plotly) {
      window.Plotly.downloadImage(plotElement, {
        format: 'png',
        width: 1200,
        height: 800,
        filename: `${selectedVariable}_plot`,
      });
    }
  };

  // Get file size in readable format
  const formatFileSize = (bytes: number): string => {
    const sizes = ['B', 'KB', 'MB', 'GB'];
    if (bytes === 0) return '0 B';
    const i = Math.floor(Math.log(bytes) / Math.log(1024));
    return `${(bytes / Math.pow(1024, i)).toFixed(2)} ${sizes[i]}`;
  };

  return (
    <Container maxWidth="xl" sx={{ mt: 3, mb: 4 }}>
      {/* Header */}
      <Box sx={{ mb: 3 }}>
        <Breadcrumbs sx={{ mb: 2 }}>
          <Link
            component="button"
            variant="body1"
            onClick={() => navigate('/jobs')}
            sx={{ cursor: 'pointer' }}
          >
            Jobs
          </Link>
          {jobId && (
            <Link
              component="button"
              variant="body1"
              onClick={() => navigate(`/jobs/${jobId}`)}
              sx={{ cursor: 'pointer' }}
            >
              {jobId.substring(0, 8)}
            </Link>
          )}
          <Typography color="text.primary">Visualization</Typography>
        </Breadcrumbs>

        <Typography variant="h4" gutterBottom>
          <InsertChartIcon sx={{ mr: 1, verticalAlign: 'bottom' }} />
          Data Visualization
        </Typography>
        <Typography variant="body1" color="text.secondary">
          Explore and visualize FIDASIM simulation output data interactively
        </Typography>
      </Box>

      <Grid container spacing={3}>
        {/* File Selection */}
        <Grid item xs={12}>
          <Paper sx={{ p: 3 }}>
            <Box display="flex" alignItems="center" justifyContent="space-between" mb={2}>
              <Typography variant="h6">
                <FolderOpenIcon sx={{ mr: 1, verticalAlign: 'bottom' }} />
                Select Data File
              </Typography>
              <Box>
                <Tooltip title="Clear cache">
                  <IconButton onClick={handleClearCache} size="small">
                    <RefreshIcon />
                  </IconButton>
                </Tooltip>
              </Box>
            </Box>

            {filesLoading ? (
              <Box display="flex" justifyContent="center" p={3}>
                <CircularProgress />
              </Box>
            ) : filesError ? (
              <Alert severity="error">Failed to load files</Alert>
            ) : files && files.length > 0 ? (
              <FormControl fullWidth>
                <InputLabel>Output File</InputLabel>
                <Select
                  value={selectedFile}
                  onChange={handleFileSelect}
                  label="Output File"
                >
                  {files.map((file) => (
                    <MenuItem key={file.path} value={file.path}>
                      <Box display="flex" alignItems="center" justifyContent="space-between" width="100%">
                        <span>{file.name}</span>
                        <Chip
                          label={formatFileSize(file.size)}
                          size="small"
                          variant="outlined"
                        />
                      </Box>
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            ) : (
              <Alert severity="info">
                {jobId
                  ? 'No output files available for this job'
                  : 'Select a job to view its output files'}
              </Alert>
            )}

            {datasetLoading && (
              <Box display="flex" justifyContent="center" p={3}>
                <CircularProgress size={24} />
                <Typography sx={{ ml: 2 }}>Loading dataset information...</Typography>
              </Box>
            )}

            {datasetError && (
              <Alert severity="error" sx={{ mt: 2 }}>
                Failed to load dataset information
              </Alert>
            )}

            {datasetInfo && (
              <Box sx={{ mt: 2 }}>
                <Typography variant="body2" color="text.secondary">
                  Dataset contains {Object.keys(datasetInfo.data_variables).length} variables
                  and {Object.keys(datasetInfo.dimensions).length} dimensions
                </Typography>
              </Box>
            )}
          </Paper>
        </Grid>

        {/* Coordinate Selection */}
        {selectedFile && (
          <Grid item xs={12} lg={4}>
            <ErrorBoundary
              onReset={() => {
                setSelectedVariable('');
                setXDimension('');
                setYDimension(null);
                setFixedCoords({});
                setIntegrateDims([]);
              }}
            >
              <CoordinateSelector
                datasetInfo={datasetInfo || null}
                selectedVariable={selectedVariable}
                onVariableChange={setSelectedVariable}
                onDimensionsChange={handleDimensionsChange}
                onFixedCoordsChange={handleFixedCoordsChange}
                onIntegrateDimsChange={handleIntegrateDimsChange}
                onVisualize={handleVisualize}
              />
            </ErrorBoundary>
          </Grid>
        )}

        {/* Plot Display */}
        {selectedFile && (
          <Grid item xs={12} lg={selectedFile ? 8 : 12}>
            <Paper sx={{ p: 2 }}>
              <Box display="flex" alignItems="center" justifyContent="space-between" mb={2}>
                <Typography variant="h6">
                  Interactive Plot
                </Typography>
                {plotData && (
                  <Tooltip title="Download plot as PNG">
                    <IconButton onClick={handleDownloadPlot} size="small">
                      <DownloadIcon />
                    </IconButton>
                  </Tooltip>
                )}
              </Box>
              <ErrorBoundary
                onReset={() => {
                  setPlotData(null);
                  setPlotError('');
                }}
              >
                <InteractivePlot
                  data={plotData}
                  loading={plotLoading}
                  error={plotError}
                  height={600}
                />
              </ErrorBoundary>
            </Paper>
          </Grid>
        )}
      </Grid>

      {/* Instructions when no job is selected */}
      {!jobId && (
        <Paper sx={{ p: 4, mt: 3, textAlign: 'center' }}>
          <InsertChartIcon sx={{ fontSize: 64, color: 'text.secondary', mb: 2 }} />
          <Typography variant="h6" gutterBottom>
            No Job Selected
          </Typography>
          <Typography color="text.secondary" paragraph>
            Navigate to a job from the Jobs page to visualize its output data
          </Typography>
          <Button
            variant="contained"
            onClick={() => navigate('/jobs')}
            sx={{ mt: 2 }}
          >
            Go to Jobs
          </Button>
        </Paper>
      )}
    </Container>
  );
};

// Add Plotly types to window
declare global {
  interface Window {
    Plotly: any;
  }
}

export default Visualization;