/**
 * Component for selecting coordinates and dimensions for visualization
 */

import React, { useState, useEffect } from 'react';
import {
  Box,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  TextField,
  Typography,
  Paper,
  Grid,
  Slider,
  Chip,
  Button,
  SelectChangeEvent,
  Tooltip,
  Switch,
  FormControlLabel,
  Checkbox,
} from '@mui/material';
import { DatasetInfo } from '../services/visualization';
import InfoIcon from '@mui/icons-material/Info';
import SwapHorizIcon from '@mui/icons-material/SwapHoriz';

interface CoordinateSelectorProps {
  datasetInfo: DatasetInfo | null;
  selectedVariable: string;
  onVariableChange: (variable: string) => void;
  onDimensionsChange: (xDim: string, yDim: string | null) => void;
  onFixedCoordsChange: (coords: Record<string, number>) => void;
  onIntegrateDimsChange: (dims: string[]) => void;
  onVisualize: () => void;
}

const CoordinateSelector: React.FC<CoordinateSelectorProps> = ({
  datasetInfo,
  selectedVariable,
  onVariableChange,
  onDimensionsChange,
  onFixedCoordsChange,
  onIntegrateDimsChange,
  onVisualize,
}) => {
  const [xDimension, setXDimension] = useState<string>('');
  const [yDimension, setYDimension] = useState<string>('');
  const [is2D, setIs2D] = useState<boolean>(false);
  const [fixedCoords, setFixedCoords] = useState<Record<string, number>>({});
  const [integrateDims, setIntegrateDims] = useState<string[]>([]);
  const [availableDims, setAvailableDims] = useState<string[]>([]);

  // Reset state when dataset changes (file switch)
  useEffect(() => {
    setXDimension('');
    setYDimension('');
    setIs2D(false);
    setFixedCoords({});
    setIntegrateDims([]);
    setAvailableDims([]);
  }, [datasetInfo]);

  useEffect(() => {
    if (datasetInfo && selectedVariable && datasetInfo.data_variables?.[selectedVariable]) {
      const varInfo = datasetInfo.data_variables[selectedVariable];
      if (varInfo?.dims && Array.isArray(varInfo.dims)) {
        setAvailableDims(varInfo.dims);

        // Initialize fixed coordinates for all dimensions
        const initialFixed: Record<string, number> = {};
        varInfo.dims.forEach((dim) => {
          if (dim !== xDimension && dim !== yDimension) {
            initialFixed[dim] = 0;
          }
        });
        setFixedCoords(initialFixed);
      }
    } else {
      // Clear state if variable is not valid
      setAvailableDims([]);
      setFixedCoords({});
    }
  }, [datasetInfo, selectedVariable, xDimension, yDimension]);

  const handleVariableChange = (event: SelectChangeEvent<string>) => {
    const variable = event.target.value;
    onVariableChange(variable);

    // Reset dimension selections
    setXDimension('');
    setYDimension('');
    setIs2D(false);
  };

  const handleXDimensionChange = (event: SelectChangeEvent<string>) => {
    const dim = event.target.value;
    setXDimension(dim);
    onDimensionsChange(dim, is2D ? yDimension : null);
  };

  const handleYDimensionChange = (event: SelectChangeEvent<string>) => {
    const dim = event.target.value;
    setYDimension(dim);
    onDimensionsChange(xDimension, dim);
  };

  const handleFixedCoordChange = (dim: string, value: number) => {
    const newFixed = { ...fixedCoords, [dim]: value };
    delete newFixed[xDimension];
    delete newFixed[yDimension];
    setFixedCoords(newFixed);
    onFixedCoordsChange(newFixed);
  };

  const handleIntegrateToggle = (dim: string) => {
    const newIntegrateDims = integrateDims.includes(dim)
      ? integrateDims.filter(d => d !== dim)
      : [...integrateDims, dim];
    setIntegrateDims(newIntegrateDims);
    onIntegrateDimsChange(newIntegrateDims);
  };

  const handleSwapDimensions = () => {
    const temp = xDimension;
    setXDimension(yDimension);
    setYDimension(temp);
    onDimensionsChange(yDimension, temp);
  };

  const handle2DToggle = (event: React.ChangeEvent<HTMLInputElement>) => {
    setIs2D(event.target.checked);
    if (!event.target.checked) {
      setYDimension('');
      onDimensionsChange(xDimension, null);
    }
  };

  if (!datasetInfo) {
    return (
      <Paper sx={{ p: 3 }}>
        <Typography color="text.secondary">
          No dataset loaded. Select a file to visualize.
        </Typography>
      </Paper>
    );
  }

  const variables = datasetInfo?.data_variables ? Object.keys(datasetInfo.data_variables) : [];
  const currentVarInfo = selectedVariable && datasetInfo?.data_variables?.[selectedVariable]
    ? datasetInfo.data_variables[selectedVariable]
    : null;

  // Get dimensions that need fixing (not selected as x or y)
  const dimsToFix = availableDims.filter(
    (dim) => dim !== xDimension && (!is2D || dim !== yDimension)
  );

  return (
    <Paper sx={{ p: 3 }}>
      <Typography variant="h6" gutterBottom>
        Visualization Settings
      </Typography>

      <Grid container spacing={3}>
        {/* Variable Selection */}
        <Grid item xs={12}>
          <FormControl fullWidth>
            <InputLabel>Variable</InputLabel>
            <Select
              value={selectedVariable}
              onChange={handleVariableChange}
              label="Variable"
            >
              {variables.map((variable) => (
                <MenuItem key={variable} value={variable}>
                  <Box display="flex" alignItems="center" justifyContent="space-between" width="100%">
                    <span>{variable}</span>
                    <Box display="flex" gap={1}>
                      {datasetInfo.data_variables?.[variable]?.dims?.map((dim) => (
                        <Chip
                          key={dim}
                          label={`${dim}: ${datasetInfo.dimensions?.[dim] || 'N/A'}`}
                          size="small"
                          variant="outlined"
                        />
                      ))}
                    </Box>
                  </Box>
                </MenuItem>
              ))}
            </Select>
          </FormControl>
          {currentVarInfo && currentVarInfo.attrs.description && (
            <Typography variant="caption" color="text.secondary" sx={{ mt: 1 }}>
              {currentVarInfo.attrs.description}
            </Typography>
          )}
        </Grid>

        {/* Plot Type Toggle */}
        {selectedVariable && availableDims.length >= 2 && (
          <Grid item xs={12}>
            <FormControlLabel
              control={
                <Switch
                  checked={is2D}
                  onChange={handle2DToggle}
                />
              }
              label="2D Plot (Heatmap)"
            />
          </Grid>
        )}

        {/* Dimension Selection */}
        {selectedVariable && (
          <>
            <Grid item xs={is2D ? 5 : 12}>
              <FormControl fullWidth>
                <InputLabel>X Axis</InputLabel>
                <Select
                  value={xDimension}
                  onChange={handleXDimensionChange}
                  label="X Axis"
                >
                  {availableDims.map((dim) => (
                    <MenuItem key={dim} value={dim} disabled={dim === yDimension}>
                      <Box display="flex" alignItems="center" justifyContent="space-between" width="100%">
                        <span>{dim}</span>
                        <Chip
                          label={`Size: ${datasetInfo.dimensions?.[dim] || 0}`}
                          size="small"
                          color="primary"
                          variant="outlined"
                        />
                      </Box>
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </Grid>

            {is2D && (
              <>
                <Grid item xs={2} display="flex" alignItems="center" justifyContent="center">
                  <Tooltip title="Swap dimensions">
                    <Button
                      onClick={handleSwapDimensions}
                      disabled={!xDimension || !yDimension}
                      variant="outlined"
                      size="small"
                    >
                      <SwapHorizIcon />
                    </Button>
                  </Tooltip>
                </Grid>

                <Grid item xs={5}>
                  <FormControl fullWidth>
                    <InputLabel>Y Axis</InputLabel>
                    <Select
                      value={yDimension}
                      onChange={handleYDimensionChange}
                      label="Y Axis"
                    >
                      {availableDims.map((dim) => (
                        <MenuItem key={dim} value={dim} disabled={dim === xDimension}>
                          <Box display="flex" alignItems="center" justifyContent="space-between" width="100%">
                            <span>{dim}</span>
                            <Chip
                              label={`Size: ${datasetInfo.dimensions[dim]}`}
                              size="small"
                              color="primary"
                              variant="outlined"
                            />
                          </Box>
                        </MenuItem>
                      ))}
                    </Select>
                  </FormControl>
                </Grid>
              </>
            )}
          </>
        )}

        {/* Fixed Coordinates */}
        {selectedVariable && dimsToFix.length > 0 && (
          <Grid item xs={12}>
            <Typography variant="subtitle1" gutterBottom>
              Fixed Coordinates
              <Tooltip title="Select indices for dimensions not used as axes">
                <InfoIcon fontSize="small" sx={{ ml: 1, color: 'text.secondary' }} />
              </Tooltip>
            </Typography>
            <Box sx={{ mt: 2 }}>
              {dimsToFix.map((dim) => {
                const maxValue = (datasetInfo.dimensions?.[dim] || 1) - 1;
                const value = fixedCoords[dim] || 0;
                const isIntegrating = integrateDims.includes(dim);

                return (
                  <Box key={dim} sx={{ mb: 3 }}>
                    <Box display="flex" alignItems="center" justifyContent="space-between" mb={1}>
                      <Typography>
                        {dim}: {isIntegrating ? 'Integrate' : `${value} / ${maxValue}`}
                      </Typography>
                      <FormControlLabel
                        control={
                          <Checkbox
                            checked={isIntegrating}
                            onChange={() => handleIntegrateToggle(dim)}
                            size="small"
                          />
                        }
                        label="Integrate"
                      />
                    </Box>
                    <Grid container spacing={2} alignItems="center">
                      <Grid item xs>
                        <Slider
                          value={value}
                          onChange={(_, newValue) =>
                            handleFixedCoordChange(dim, newValue as number)
                          }
                          min={0}
                          max={maxValue}
                          step={1}
                          marks
                          valueLabelDisplay="auto"
                          disabled={isIntegrating}
                        />
                      </Grid>
                      <Grid item>
                        <TextField
                          value={value}
                          onChange={(e) =>
                            handleFixedCoordChange(dim, parseInt(e.target.value) || 0)
                          }
                          type="number"
                          inputProps={{
                            min: 0,
                            max: maxValue,
                          }}
                          size="small"
                          sx={{ width: 80 }}
                          disabled={isIntegrating}
                        />
                      </Grid>
                    </Grid>
                  </Box>
                );
              })}
            </Box>
          </Grid>
        )}

        {/* Visualize Button */}
        <Grid item xs={12}>
          <Button
            variant="contained"
            color="primary"
            onClick={onVisualize}
            disabled={!selectedVariable || !xDimension}
            fullWidth
            size="large"
          >
            Visualize Data
          </Button>
        </Grid>
      </Grid>
    </Paper>
  );
};

export default CoordinateSelector;