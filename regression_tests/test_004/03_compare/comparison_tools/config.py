"""Read Test 004 comparison settings and construct matched file pairs."""

from dataclasses import dataclass
import math
from pathlib import Path

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    read_namelist,
    require_boolean,
    require_choice,
    require_existing_file,
    require_integer,
    require_real,
    require_string,
    validate_schema,
)
from test_004_tools import discover_distributions, read_test_config


CONFIG_SCHEMA = {
    "compare": {
        "required": True,
        "required_fields": [
            "test_config",
            "rate_relative_tolerance",
            "rate_sigma_tolerance",
        ],
        "optional_fields": ["comment", "generate_plots"],
    },
    "plot_data_block": {
        "required": False,
        "required_fields": [],
        "optional_fields": [
            "scale",
            "emax",
            "fmin",
            "fmax",
            "enable_colorbar",
            "colormap",
            "contour_levels",
        ],
    },
}

SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]


@dataclass(frozen=True)
class FilePair:
    """Paths for one deterministic and Monte Carlo result pair."""

    case_index: int
    source_distribution: Path
    deterministic: Path
    monte_carlo: Path
    expected_provenance: dict

    @property
    def basename(self):
        """Return the shared HDF5 basename."""
        return self.deterministic.name

    @property
    def stem(self):
        """Return the shared filename stem."""
        return self.deterministic.stem


@dataclass
class ComparisonInputs:
    """Resolved upstream inputs needed by the comparison workflow."""

    file_pairs: list
    output_directory: Path
    test_config: dict


def _repository_relative(path):
    """Return repository-relative provenance when the path is inside the repo."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _optional_comment(block):
    """Return the optional human-readable comparison comment."""
    if "comment" not in block:
        return ""
    return require_string(block["comment"], "compare comment")


def _positive_real(value, label):
    """Return a finite, strictly positive floating-point value."""
    result = require_real(value, label)
    if not math.isfinite(result) or result <= 0.0:
        raise ConfigError(f"{label} must be finite and greater than zero.")
    return result


def _normalize_plot_limit(value, label):
    """Normalize an automatic or finite numerical plot limit."""
    if value is None:
        return None
    if isinstance(value, str):
        normalized = normalize_string(value)
        if normalized == "auto":
            return normalized
        raise ConfigError(f"{label} must be a real number or 'auto'.")

    limit = require_real(value, label)
    if not math.isfinite(limit):
        raise ConfigError(f"{label} must be finite.")
    return limit


def _read_plot_config(block):
    """Return canonical comparison-plot settings."""
    scale = require_choice(
        normalize_string(block.get("scale", "lin")),
        SUPPORTED_SCALES,
        "scale",
    )
    colormap = require_choice(
        normalize_string(block.get("colormap", "viridis")),
        SUPPORTED_COLORMAPS,
        "colormap",
    )
    enable_colorbar = require_boolean(
        block.get("enable_colorbar", True),
        "enable_colorbar",
    )
    contour_levels = require_integer(
        block.get("contour_levels", 30),
        "contour_levels",
    )
    if contour_levels < 2:
        raise ConfigError("contour_levels must be at least 2.")

    emax = _normalize_plot_limit(block.get("emax"), "emax")
    fmin = _normalize_plot_limit(block.get("fmin"), "fmin")
    fmax = _normalize_plot_limit(block.get("fmax"), "fmax")

    if isinstance(emax, float) and emax <= 0.0:
        raise ConfigError("emax must be greater than zero.")
    if isinstance(fmin, float) and isinstance(fmax, float) and fmin >= fmax:
        raise ConfigError("fmin must be less than fmax.")
    if scale == "log":
        if isinstance(fmin, float) and fmin <= 0.0:
            raise ConfigError("fmin must be greater than zero for log scale.")
        if isinstance(fmax, float) and fmax <= 0.0:
            raise ConfigError("fmax must be greater than zero for log scale.")

    return {
        "scale": scale,
        "emax": emax,
        "fmin": fmin,
        "fmax": fmax,
        "enable_colorbar": enable_colorbar,
        "colormap": colormap,
        "contour_levels": contour_levels,
    }


def read_comparison_config(config_filename):
    """Read and validate one Test 004 Stage 3 namelist."""
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    compare_block = blocks["compare"]
    test_config = normalize_path(
        value=require_string(compare_block["test_config"], "test_config"),
        config_path=config_path,
        field_label="test_config",
    )
    require_existing_file(test_config, "test_config")

    generate_plots = require_boolean(
        compare_block.get("generate_plots", True),
        "generate_plots",
    )
    plot_block = blocks.get("plot_data_block", {})

    return {
        "compare": {
            "comment": _optional_comment(compare_block),
            "config_path": config_path,
            "test_config": test_config,
            "rate_relative_tolerance": _positive_real(
                compare_block["rate_relative_tolerance"],
                "rate_relative_tolerance",
            ),
            "rate_sigma_tolerance": _positive_real(
                compare_block["rate_sigma_tolerance"],
                "rate_sigma_tolerance",
            ),
            "generate_plots": generate_plots,
        },
        "plot_data_block": _read_plot_config(plot_block),
    }


def build_file_pairs(test_config_filename):
    """Derive and preflight all files and settings consumed by Stage 3.

    The unified Test 004 configuration identifies the smooth input
    distributions and the common output root. Each source distribution stem
    maps to one deterministic and one Monte Carlo ion-sink filename.
    """
    test_config = read_test_config(test_config_filename)
    if not test_config["deterministic"]["save_data"]:
        raise ConfigError(
            "The selected Test 004 configuration disables deterministic "
            "HDF5 output."
        )
    if not test_config["monte_carlo"]["save_data"]:
        raise ConfigError(
            "The selected Test 004 configuration disables Monte Carlo "
            "HDF5 output."
        )

    provenance = discover_distributions(
        test_config["test_case"]["input_distribution_config"]
    )
    deterministic_directory = test_config["save_data_block"][
        "deterministic_directory"
    ]
    monte_carlo_directory = test_config["save_data_block"][
        "monte_carlo_directory"
    ]
    common_provenance = {
        "input_distribution_config": _repository_relative(
            provenance["run_config"]
        ),
        "atomic_tables_file": _repository_relative(
            test_config["test_case"]["tables_filename"]
        ),
        "test_case_config": _repository_relative(
            test_config["test_case"]["config_path"]
        ),
        "test_case_comment": test_config["test_case"]["comment"],
        "level_split_method": test_config["neutrals"][
            "level_split_method"
        ],
    }

    pairs = []
    for case in provenance["cases"]:
        stem = f"{Path(case['input_path']).stem}_ion_sink"
        expected_provenance = {
            **common_provenance,
            "source_distribution": _repository_relative(
                case["input_path"]
            ),
        }
        pairs.append(
            FilePair(
                case_index=case["index"],
                source_distribution=Path(case["input_path"]),
                deterministic=deterministic_directory / f"{stem}.h5",
                monte_carlo=monte_carlo_directory / f"{stem}.h5",
                expected_provenance=expected_provenance,
            )
        )

    # Check the complete collection before reading or plotting any pair. This
    # avoids producing a partial comparison when one upstream stage is absent.
    missing_files = []
    for pair in pairs:
        if not pair.deterministic.is_file():
            missing_files.append(str(pair.deterministic))
        if not pair.monte_carlo.is_file():
            missing_files.append(str(pair.monte_carlo))
    if missing_files:
        raise ConfigError(
            "Missing Test 004 result files: " + ", ".join(missing_files)
        )

    output_directory = (
        test_config["save_data_block"]["output_directory"] / "comparison"
    )
    return ComparisonInputs(
        file_pairs=pairs,
        output_directory=output_directory,
        test_config=test_config,
    )
