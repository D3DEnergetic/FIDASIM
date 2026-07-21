"""Select readers that produce the canonical energy-pitch representation.

Every registered reader returns ``z, r, pitch, energy, f, denf``, with ``f``
ordered as ``(z, r, pitch, energy)`` and ``denf`` ordered as ``(z, r)``.
"""

from pathlib import Path

from .fidasim_h5 import load_fidasim_h5_distribution


def load_input_distribution(config):
    """Select a reader and load the canonical distribution and density.

    Args:
        config (dict): Validated Stage 1 configuration organized by namelist
        block.

    Returns:
        tuple: ``z, r, pitch, energy, f, denf`` from the selected reader.
        ``f`` is ordered as ``(z, r, pitch, energy)`` and ``denf`` as
        ``(z, r)``.

    Raises:
        FileNotFoundError: If the configured input file does not exist.
        NotImplementedError: If the recognized input type has no reader.
        ValueError: If no reader is available for the configured input type.
    """
    # Extract the input file information from the validated configuration.
    input_config = config["input"]
    input_file_type = input_config["input_file_type"]
    input_filename = input_config["input_filename"]
    input_path = Path(input_filename)

    # Select the format-specific reader. Add future input types as new branches.
    if input_file_type == "fidasim_h5":
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        return load_fidasim_h5_distribution(input_path)
    elif input_file_type == "cql3d_f4d":
        raise NotImplementedError(
            "The input distribution reader for 'cql3d_f4d' has not been "
            "implemented."
        )

    raise ValueError(f"No reader is available for input type '{input_file_type}'.")
