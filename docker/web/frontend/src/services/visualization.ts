/**
 * Visualization service for interacting with the backend API
 */

import axios from 'axios';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000';

export interface DatasetInfo {
  dimensions: Record<string, number>;
  coordinates: Record<string, {
    dims: string[];
    shape: number[];
    dtype: string;
    attrs: Record<string, any>;
  }>;
  data_variables: Record<string, {
    dims: string[];
    shape: number[];
    dtype: string;
    attrs: Record<string, any>;
  }>;
  attributes: Record<string, any>;
}

export interface VariableData {
  variable: string;
  dims: string[];
  shape: number[];
  data: any;
  coords: Record<string, any>;
  attrs: Record<string, any>;
}

export interface PlotData1D {
  type: '1D';
  x: number[];
  y: number[];
  x_label: string;
  y_label: string;
  x_units: string;
  y_units: string;
}

export interface PlotData2D {
  type: '2D';
  x: number[];
  y: number[];
  z: number[][];
  x_label: string;
  y_label: string;
  z_label: string;
  x_units: string;
  y_units: string;
  z_units: string;
}

export type PlotData = PlotData1D | PlotData2D;

export interface VisualizationFile {
  name: string;
  path: string;
  size: number;
  modified: number;
}

class VisualizationService {
  private getAuthHeaders() {
    const token = localStorage.getItem('token');
    return {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    };
  }

  async listFiles(jobId: string): Promise<VisualizationFile[]> {
    const response = await axios.get(
      `${API_URL}/api/visualization/files/${jobId}`,
      this.getAuthHeaders()
    );
    return response.data;
  }

  async getDatasetInfo(filePath: string): Promise<DatasetInfo> {
    const response = await axios.get(
      `${API_URL}/api/visualization/dataset/info`,
      {
        ...this.getAuthHeaders(),
        params: { file_path: filePath },
      }
    );
    return response.data;
  }

  async listVariables(filePath: string): Promise<{
    data_variables: string[];
    coordinates: string[];
    dimensions: string[];
  }> {
    const response = await axios.get(
      `${API_URL}/api/visualization/dataset/variables`,
      {
        ...this.getAuthHeaders(),
        params: { file_path: filePath },
      }
    );
    return response.data;
  }

  async getVariableData(
    filePath: string,
    variable: string,
    slices?: Record<string, any>
  ): Promise<VariableData> {
    const response = await axios.post(
      `${API_URL}/api/visualization/dataset/variable`,
      {
        variable,
        slices,
      },
      {
        ...this.getAuthHeaders(),
        params: { file_path: filePath },
      }
    );
    return response.data;
  }

  async getPlotData(
    filePath: string,
    variable: string,
    xDim: string,
    yDim?: string,
    fixedCoords?: Record<string, number>,
    integrateDims?: string[]
  ): Promise<PlotData> {
    const response = await axios.post(
      `${API_URL}/api/visualization/plot/data`,
      {
        variable,
        x_dim: xDim,
        y_dim: yDim,
        fixed_coords: fixedCoords,
        integrate_dims: integrateDims,
      },
      {
        ...this.getAuthHeaders(),
        params: { file_path: filePath },
      }
    );
    return response.data;
  }

  async clearCache(filePath?: string): Promise<void> {
    await axios.post(
      `${API_URL}/api/visualization/cache/clear`,
      {},
      {
        ...this.getAuthHeaders(),
        params: filePath ? { file_path: filePath } : {},
      }
    );
  }
}

export default new VisualizationService();