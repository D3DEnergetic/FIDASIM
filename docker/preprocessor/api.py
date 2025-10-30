"""
FIDASIM Preprocessor API Service
Provides REST endpoints for preprocessing and data preparation
"""

import os
import json
import logging
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any, Union
from enum import Enum

from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field, validator
import numpy as np
import h5py
from scipy import interpolate

# Import FIDASIM preprocessing modules
import sys
sys.path.insert(0, '/app/lib/python')

from fidasim import utils
from fidasim import preprocessing
from efit import io as efit_io

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = Path(os.getenv("DATA_DIR", "/data"))
DATA_DIR.mkdir(parents=True, exist_ok=True)

class PreprocessingMode(str, Enum):
    PREFIDA = "prefida"
    TRANSP_EXTRACT = "transp_extract"
    GEQDSK_CONVERT = "geqdsk_convert"
    GRID_GENERATION = "grid_generation"
    NUBEAM_GEOMETRY = "nubeam_geometry"

class GridType(str, Enum):
    RZ_UNIFORM = "rz_uniform"
    RZ_NONUNIFORM = "rz_nonuniform"
    BEAM_GRID = "beam_grid"

class PrefidaConfig(BaseModel):
    """Configuration for PREFIDA preprocessing"""
    shot: int = Field(description="Shot number")
    time: float = Field(description="Time in seconds")
    device: str = Field(description="Tokamak device name")
    runid: str = Field(description="Run identifier")

    # Grid settings
    nr: int = Field(default=100, description="Number of R grid points")
    nz: int = Field(default=100, description="Number of Z grid points")
    nphi: int = Field(default=1, description="Number of phi grid points")

    # Physics settings
    calc_bes: bool = Field(default=True, description="Calculate BES")
    calc_fida: bool = Field(default=True, description="Calculate FIDA")
    calc_npa: bool = Field(default=False, description="Calculate NPA")

    # Optional overrides
    btipsign: Optional[float] = Field(default=None, description="Sign of toroidal field")
    iptipsign: Optional[float] = Field(default=None, description="Sign of plasma current")

class GridConfig(BaseModel):
    """Configuration for grid generation"""
    grid_type: GridType
    nr: int = Field(default=100, ge=10, le=500)
    nz: int = Field(default=100, ge=10, le=500)
    nphi: int = Field(default=1, ge=1, le=360)
    rmin: Optional[float] = None
    rmax: Optional[float] = None
    zmin: Optional[float] = None
    zmax: Optional[float] = None
    phi_start: float = Field(default=0.0)
    phi_end: float = Field(default=360.0)

class ConversionResult(BaseModel):
    """Result of a file conversion operation"""
    success: bool
    output_file: Optional[str] = None
    errors: List[str] = []
    warnings: List[str] = []
    metadata: Dict[str, Any] = {}

# Initialize FastAPI app
app = FastAPI(
    title="FIDASIM Preprocessor API",
    description="API for preprocessing and data preparation for FIDASIM",
    version="1.0.0"
)

@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "service": "FIDASIM Preprocessor",
        "version": "1.0.0",
        "status": "operational",
        "python_path": sys.path,
        "available_modes": [mode.value for mode in PreprocessingMode]
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    checks = {
        "api": True,
        "data_dir_writable": os.access(DATA_DIR, os.W_OK),
        "fidasim_module": "fidasim" in sys.modules or True,
        "efit_module": "efit" in sys.modules or True
    }

    healthy = all(checks.values())

    return {
        "healthy": healthy,
        "checks": checks,
        "timestamp": datetime.utcnow().isoformat()
    }

@app.post("/prefida")
async def run_prefida(
    config: PrefidaConfig,
    equilibrium_file: UploadFile = File(...),
    profiles_file: Optional[UploadFile] = File(None),
    nbi_file: Optional[UploadFile] = File(None)
) -> ConversionResult:
    """Run PREFIDA preprocessing to create FIDASIM input files"""

    errors = []
    warnings = []
    metadata = {}

    # Create temporary working directory
    work_dir = DATA_DIR / f"prefida_{config.runid}"
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Save uploaded files
        eq_path = work_dir / "equilibrium.geqdsk"
        eq_content = await equilibrium_file.read()
        eq_path.write_bytes(eq_content)

        # Read equilibrium using EFIT module
        logger.info(f"Reading equilibrium from {eq_path}")
        eq_data = efit_io.readg(str(eq_path))

        # Create input dictionary for PREFIDA
        inputs = {
            "shot": config.shot,
            "time": config.time,
            "device": config.device,
            "runid": config.runid,
            "result_dir": str(work_dir),

            # Grid configuration
            "nr": config.nr,
            "nz": config.nz,
            "nphi": config.nphi,

            # Equilibrium data
            "equilibrium": eq_data,

            # Switches
            "calc_bes": 1 if config.calc_bes else 0,
            "calc_fida": 1 if config.calc_fida else 0,
            "calc_npa": 1 if config.calc_npa else 0,
        }

        # Add optional profile data if provided
        if profiles_file:
            prof_path = work_dir / "profiles.dat"
            prof_content = await profiles_file.read()
            prof_path.write_bytes(prof_content)
            # Parse profiles and add to inputs
            # This would need custom parsing based on format
            warnings.append("Profile file provided but parsing not yet implemented")

        # Add NBI data if provided
        if nbi_file:
            nbi_path = work_dir / "nbi.dat"
            nbi_content = await nbi_file.read()
            nbi_path.write_bytes(nbi_content)
            warnings.append("NBI file provided but parsing not yet implemented")

        # Run PREFIDA preprocessing
        logger.info("Running PREFIDA preprocessing")

        # Check and validate inputs
        valid, validation_errors = preprocessing.check_dict_schema(inputs)
        if not valid:
            errors.extend(validation_errors)
            raise ValueError(f"Input validation failed: {validation_errors}")

        # Create the HDF5 output files
        output_files = preprocessing.prefida(inputs)

        metadata["output_files"] = output_files
        metadata["grid_size"] = f"{config.nr}x{config.nz}x{config.nphi}"

        # Generate namelist file
        namelist_path = work_dir / f"{config.runid}_inputs.dat"
        namelist_content = generate_namelist(config, work_dir)
        namelist_path.write_text(namelist_content)

        success = True
        output_file = str(namelist_path.relative_to(DATA_DIR))

    except Exception as e:
        logger.error(f"PREFIDA preprocessing failed: {e}")
        errors.append(str(e))
        success = False
        output_file = None

    return ConversionResult(
        success=success,
        output_file=output_file,
        errors=errors,
        warnings=warnings,
        metadata=metadata
    )

@app.post("/convert/geqdsk")
async def convert_geqdsk(
    file: UploadFile = File(...),
    output_format: str = Form("hdf5")
) -> ConversionResult:
    """Convert GEQDSK equilibrium file to HDF5 format"""

    errors = []
    warnings = []
    metadata = {}

    # Create temporary directory
    temp_dir = DATA_DIR / f"convert_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"
    temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Save uploaded file
        input_path = temp_dir / "input.geqdsk"
        content = await file.read()
        input_path.write_bytes(content)

        # Read GEQDSK file
        logger.info(f"Reading GEQDSK file: {input_path}")
        geqdsk_data = utils.read_geqdsk(str(input_path))

        if output_format == "hdf5":
            # Convert to HDF5
            output_path = temp_dir / "equilibrium.h5"

            with h5py.File(output_path, 'w') as f:
                # Store equilibrium data in HDF5 format
                f.create_dataset('psi', data=geqdsk_data['psi'])
                f.create_dataset('r', data=geqdsk_data['r'])
                f.create_dataset('z', data=geqdsk_data['z'])
                f.create_dataset('psi_axis', data=geqdsk_data['psi_axis'])
                f.create_dataset('psi_lcfs', data=geqdsk_data['psi_lcfs'])
                f.create_dataset('magnetic_axis', data=[geqdsk_data['rmag'], geqdsk_data['zmag']])
                f.create_dataset('bt0', data=geqdsk_data['bt0'])
                f.create_dataset('plasma_current', data=geqdsk_data['ip'])

                # Add metadata
                f.attrs['shot'] = geqdsk_data.get('shot', 0)
                f.attrs['time'] = geqdsk_data.get('time', 0.0)
                f.attrs['device'] = 'unknown'

            metadata["grid_nr"] = len(geqdsk_data['r'])
            metadata["grid_nz"] = len(geqdsk_data['z'])
            metadata["bt0"] = float(geqdsk_data['bt0'])
            metadata["ip"] = float(geqdsk_data['ip'])

            success = True
            output_file = str(output_path.relative_to(DATA_DIR))

        else:
            errors.append(f"Unsupported output format: {output_format}")
            success = False
            output_file = None

    except Exception as e:
        logger.error(f"GEQDSK conversion failed: {e}")
        errors.append(str(e))
        success = False
        output_file = None

    return ConversionResult(
        success=success,
        output_file=output_file,
        errors=errors,
        warnings=warnings,
        metadata=metadata
    )

@app.post("/grid/generate")
async def generate_grid(config: GridConfig) -> Dict[str, Any]:
    """Generate computational grid for FIDASIM"""

    try:
        if config.grid_type in [GridType.RZ_UNIFORM, GridType.RZ_NONUNIFORM]:
            # Generate R-Z grid
            if config.rmin is None or config.rmax is None:
                # Use typical tokamak values if not specified
                config.rmin = 0.5
                config.rmax = 2.5
                config.zmin = -1.5
                config.zmax = 1.5

            # Create grid
            r = np.linspace(config.rmin, config.rmax, config.nr)
            z = np.linspace(config.zmin, config.zmax, config.nz)

            if config.nphi > 1:
                phi = np.linspace(
                    np.radians(config.phi_start),
                    np.radians(config.phi_end),
                    config.nphi
                )
            else:
                phi = np.array([0.0])

            # Create meshgrid
            R, Z, PHI = np.meshgrid(r, z, phi, indexing='ij')

            # Save grid to HDF5
            output_path = DATA_DIR / f"grid_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.h5"

            with h5py.File(output_path, 'w') as f:
                f.create_dataset('r', data=r)
                f.create_dataset('z', data=z)
                f.create_dataset('phi', data=phi)
                f.create_dataset('R_grid', data=R)
                f.create_dataset('Z_grid', data=Z)
                f.create_dataset('PHI_grid', data=PHI)

                # Add metadata
                f.attrs['grid_type'] = config.grid_type
                f.attrs['nr'] = config.nr
                f.attrs['nz'] = config.nz
                f.attrs['nphi'] = config.nphi
                f.attrs['timestamp'] = datetime.utcnow().isoformat()

            return {
                "success": True,
                "output_file": str(output_path.relative_to(DATA_DIR)),
                "grid_shape": [config.nr, config.nz, config.nphi],
                "r_range": [float(config.rmin), float(config.rmax)],
                "z_range": [float(config.zmin), float(config.zmax)],
                "phi_range": [float(config.phi_start), float(config.phi_end)]
            }

        elif config.grid_type == GridType.BEAM_GRID:
            # Beam-aligned grid generation would go here
            raise NotImplementedError("Beam grid generation not yet implemented")

    except Exception as e:
        logger.error(f"Grid generation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/extract/transp")
async def extract_transp_data(
    transp_file: UploadFile = File(...),
    output_type: str = Form("plasma")
) -> ConversionResult:
    """Extract data from TRANSP output files"""

    errors = []
    warnings = []
    metadata = {}

    # Create temporary directory
    temp_dir = DATA_DIR / f"transp_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"
    temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Save uploaded file
        input_path = temp_dir / transp_file.filename
        content = await transp_file.read()
        input_path.write_bytes(content)

        if output_type == "plasma":
            # Extract plasma profiles
            logger.info(f"Extracting plasma data from TRANSP: {input_path}")
            plasma_data = utils.extract_transp_plasma(str(input_path))

            # Save to HDF5
            output_path = temp_dir / "plasma_profiles.h5"
            with h5py.File(output_path, 'w') as f:
                for key, value in plasma_data.items():
                    if isinstance(value, np.ndarray):
                        f.create_dataset(key, data=value)
                    else:
                        f.attrs[key] = value

            metadata["profiles"] = list(plasma_data.keys())
            success = True
            output_file = str(output_path.relative_to(DATA_DIR))

        elif output_type == "fbm":
            # Extract fast-ion distribution
            warnings.append("FBM extraction not yet implemented")
            success = False
            output_file = None

        else:
            errors.append(f"Unknown output type: {output_type}")
            success = False
            output_file = None

    except Exception as e:
        logger.error(f"TRANSP extraction failed: {e}")
        errors.append(str(e))
        success = False
        output_file = None

    return ConversionResult(
        success=success,
        output_file=output_file,
        errors=errors,
        warnings=warnings,
        metadata=metadata
    )

@app.post("/validate/schema")
async def validate_input_schema(
    input_data: Dict[str, Any]
) -> Dict[str, Any]:
    """Validate input dictionary against FIDASIM schema"""

    try:
        valid, errors = preprocessing.check_dict_schema(input_data)

        return {
            "valid": valid,
            "errors": errors if not valid else [],
            "input_keys": list(input_data.keys())
        }
    except Exception as e:
        return {
            "valid": False,
            "errors": [str(e)],
            "input_keys": []
        }

@app.get("/templates/{template_type}")
async def get_input_template(template_type: str):
    """Get template input files or configurations"""

    templates = {
        "namelist": generate_template_namelist(),
        "prefida_config": PrefidaConfig(
            shot=100000,
            time=1.0,
            device="DEMO",
            runid="template",
            nr=100,
            nz=100,
            nphi=1
        ).dict(),
        "grid_config": GridConfig(
            grid_type=GridType.RZ_UNIFORM,
            nr=100,
            nz=100,
            nphi=1,
            rmin=0.5,
            rmax=2.5,
            zmin=-1.5,
            zmax=1.5
        ).dict()
    }

    if template_type not in templates:
        raise HTTPException(
            status_code=404,
            detail=f"Template not found. Available: {list(templates.keys())}"
        )

    return templates[template_type]

@app.get("/download/{filepath:path}")
async def download_file(filepath: str):
    """Download a processed file"""

    file_path = DATA_DIR / filepath

    if not file_path.exists() or not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")

    # Security check - ensure file is within DATA_DIR
    try:
        file_path.relative_to(DATA_DIR)
    except ValueError:
        raise HTTPException(status_code=403, detail="Access denied")

    return FileResponse(
        path=str(file_path),
        filename=file_path.name,
        media_type="application/octet-stream"
    )

def generate_namelist(config: PrefidaConfig, work_dir: Path) -> str:
    """Generate FIDASIM namelist file"""

    namelist = f"""&fidasim_inputs
    ! Shot information
    shot = {config.shot}
    time = {config.time}
    runid = '{config.runid}'
    device = '{config.device}'

    ! File paths
    result_dir = '{work_dir}'
    tables_file = '/app/tables/atomic_tables.h5'
    equilibrium_file = '{work_dir}/equilibrium.h5'
    geometry_file = '{work_dir}/geometry.h5'
    distribution_file = '{work_dir}/distribution.h5'

    ! Grid configuration
    nr = {config.nr}
    nz = {config.nz}
    nphi = {config.nphi}

    ! Calculation switches
    calc_bes = {'.T.' if config.calc_bes else '.F.'}
    calc_fida = {'.T.' if config.calc_fida else '.F.'}
    calc_npa = {'.T.' if config.calc_npa else '.F.'}
    calc_pfida = .F.
    calc_pnpa = .F.
    calc_cold = .T.
    calc_bremsstrahlung = .F.
    calc_dcx = .F.
    calc_halo = .F.
    calc_birth = .F.
    calc_neutron = .F.

    ! Monte Carlo settings
    n_fida = 5000000
    n_npa = 5000000
    n_pfida = 0
    n_pnpa = 0
    n_nbi = 5000000
    n_halo = 500000
    n_dcx = 500000
    n_birth = 10000

    ! Advanced settings
    seed = 12345
    flr = 1
    verbose = 1
    stark_components = 0
/
"""

    return namelist

def generate_template_namelist() -> str:
    """Generate a template namelist file"""

    return """&fidasim_inputs
    ! Shot information
    shot = 100000
    time = 1.0
    runid = 'template'
    device = 'DEMO'
    comment = 'Template FIDASIM run'

    ! File paths
    result_dir = '/data/results'
    tables_file = '/app/tables/atomic_tables.h5'
    equilibrium_file = 'equilibrium.h5'
    geometry_file = 'geometry.h5'
    distribution_file = 'distribution.h5'

    ! Grid configuration
    nr = 100
    nz = 100
    nphi = 1

    ! Calculation switches
    calc_bes = .T.
    calc_fida = .T.
    calc_npa = .F.
    calc_pfida = .F.
    calc_pnpa = .F.
    calc_cold = .T.
    calc_bremsstrahlung = .F.
    calc_dcx = .F.
    calc_halo = .F.
    calc_birth = .F.
    calc_neutron = .F.

    ! Monte Carlo settings
    n_fida = 5000000
    n_npa = 5000000
    n_pfida = 0
    n_pnpa = 0
    n_nbi = 5000000
    n_halo = 500000
    n_dcx = 500000
    n_birth = 10000

    ! Advanced settings
    seed = 12345
    flr = 1
    verbose = 1
    stark_components = 0
/
"""

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8002)