/**
 * Interactive plot component using Plotly
 */

import React from 'react';
import Plot from 'react-plotly.js';
import { Box, Paper, Typography, CircularProgress } from '@mui/material';
import { PlotData } from '../services/visualization';

interface InteractivePlotProps {
  data: PlotData | null;
  loading?: boolean;
  error?: string;
  title?: string;
  height?: number;
}

const InteractivePlot: React.FC<InteractivePlotProps> = ({
  data,
  loading = false,
  error,
  title,
  height = 500,
}) => {
  if (loading) {
    return (
      <Paper sx={{ p: 3, display: 'flex', justifyContent: 'center', alignItems: 'center', height }}>
        <CircularProgress />
      </Paper>
    );
  }

  if (error) {
    return (
      <Paper sx={{ p: 3 }}>
        <Typography color="error">Error loading plot: {error}</Typography>
      </Paper>
    );
  }

  if (!data) {
    return (
      <Paper sx={{ p: 3 }}>
        <Typography color="text.secondary">
          Select a variable and dimensions to visualize the data
        </Typography>
      </Paper>
    );
  }

  // Prepare Plotly data based on plot type
  let plotlyData: Plotly.Data[];
  let layout: Partial<Plotly.Layout>;

  if (data.type === '1D') {
    // 1D line plot
    plotlyData = [{
      x: data.x,
      y: data.y,
      type: 'scatter',
      mode: 'lines',
      name: data.y_label,
      line: {
        color: '#1976d2',
        width: 2,
      },
    }];

    layout = {
      title: title || `${data.y_label} vs ${data.x_label}`,
      xaxis: {
        title: {
          text: data.x_units ? `${data.x_label} (${data.x_units})` : data.x_label,
        },
        gridcolor: '#e0e0e0',
      },
      yaxis: {
        title: {
          text: data.y_units ? `${data.y_label} (${data.y_units})` : data.y_label,
        },
        gridcolor: '#e0e0e0',
      },
      height,
      paper_bgcolor: '#fafafa',
      plot_bgcolor: '#ffffff',
      margin: {
        l: 60,
        r: 40,
        b: 60,
        t: 60,
      },
      hovermode: 'closest',
    };
  } else {
    // 2D heatmap
    plotlyData = [{
      x: data.x,
      y: data.y,
      z: data.z,
      type: 'heatmap',
      colorscale: 'Viridis',
      colorbar: {
        title: {
          text: data.z_units ? `${data.z_label} (${data.z_units})` : data.z_label,
        },
      },
    }];

    layout = {
      title: title || `${data.z_label}`,
      xaxis: {
        title: {
          text: data.x_units ? `${data.x_label} (${data.x_units})` : data.x_label,
        },
        gridcolor: '#e0e0e0',
      },
      yaxis: {
        title: {
          text: data.y_units ? `${data.y_label} (${data.y_units})` : data.y_label,
        },
        gridcolor: '#e0e0e0',
      },
      height,
      paper_bgcolor: '#fafafa',
      plot_bgcolor: '#ffffff',
      margin: {
        l: 60,
        r: 100,
        b: 60,
        t: 60,
      },
    };
  }

  const config: Partial<Plotly.Config> = {
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['select2d', 'lasso2d'],
    toImageButtonOptions: {
      format: 'png',
      filename: title || 'fidasim_plot',
      height: 800,
      width: 1200,
      scale: 1,
    },
  };

  return (
    <Paper sx={{ p: 2 }}>
      <Box sx={{ width: '100%', overflow: 'hidden' }}>
        <Plot
          data={plotlyData}
          layout={layout}
          config={config}
          style={{ width: '100%', height: '100%' }}
          useResizeHandler
        />
      </Box>
    </Paper>
  );
};

export default InteractivePlot;