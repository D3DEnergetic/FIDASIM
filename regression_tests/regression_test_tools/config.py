"""Common operations for reading and validating namelist configurations."""

from collections.abc import Mapping
from pathlib import Path

import f90nml


class ConfigError(ValueError):
    """Identify errors caused by invalid regression-test configuration.

    This specialized ``ValueError`` allows command-line programs to catch and
    report configuration problems separately from unexpected programming,
    numerical, or file-processing errors.
    """


def print_config(config):
    """Print the keys of a nested configuration dictionary as a tree.

    Args:
        config (mapping): Nested configuration dictionary whose keys will be
        printed.

    Returns:
        None.

    Raises:
        ConfigError: If ``config`` is not a mapping.
    """
    if not isinstance(config, Mapping):
        raise ConfigError("config must be a mapping.")

    print("Configuration")
    _print_dictionary_branches(
        dictionary=config,
        prefix="",
    )


def _print_dictionary_branches(dictionary, prefix):
    """Print one level of a nested dictionary tree.

    Args:
        dictionary (mapping): Current dictionary level to print.
        prefix (str): Tree characters inherited from parent levels.

    Returns:
        None.
    """
    entries = list(dictionary.items())
    number_of_entries = len(entries)

    for entry_index, entry in enumerate(entries):
        key, value = entry
        is_last_entry = entry_index == number_of_entries - 1

        if is_last_entry:
            branch = "└── "
            child_prefix = prefix + "    "
        else:
            branch = "├── "
            child_prefix = prefix + "│   "

        print(f"{prefix}{branch}{key}")

        if isinstance(value, Mapping):
            _print_dictionary_branches(
                dictionary=value,
                prefix=child_prefix,
            )


def read_namelist(config_path):
    """Read a Fortran namelist configuration file.

    Args:
        config_path (str or Path): Path to the namelist configuration file.

    Returns:
        tuple: A two-item tuple containing the resolved configuration ``Path``
        followed by the namelist blocks returned by ``f90nml.read``.

    Raises:
        ConfigError: If ``config_path`` is invalid, does not identify a file,
        or the namelist cannot be read.
    """

    # Normalize the configuration path before opening the file. Resolving the
    # path here also provides a stable base for other configured paths.
    try:
        resolved_config_path = Path(config_path).expanduser().resolve()
    except TypeError as error:
        raise ConfigError("config_path must be a string or Path.") from error

    require_existing_file(
        path=resolved_config_path,
        field_label="Configuration file",
    )

    # Read the configuration file:
    try:
        blocks = f90nml.read(resolved_config_path)
    except (OSError, ValueError) as error:
        raise ConfigError(
            f"Unable to read namelist configuration "
            f"'{resolved_config_path}': {error}"
        ) from error

    return resolved_config_path, blocks


def require_blocks(blocks, required_blocks):
    """Check that all required namelist blocks are present.

    Args:
        blocks (mapping): Namelist blocks returned by ``f90nml.read``.
        required_blocks (list or tuple of str): Names of blocks that must be
        present. Names do not include the leading ``&``.

    Returns:
        None.

    Raises:
        ConfigError: If one or more required blocks are missing.
    """
    missing_blocks = []

    for block_name in required_blocks:
        if block_name not in blocks:
            missing_blocks.append(block_name)

    if missing_blocks:
        formatted_names = []
        for block_name in missing_blocks:
            formatted_names.append(f"&{block_name}")

        missing_text = ", ".join(formatted_names)
        raise ConfigError(f"Missing required namelist block(s): {missing_text}")


def reject_unknown_blocks(blocks, allowed_blocks):
    """Reject namelist blocks that are not part of an application schema.

    Args:
        blocks (mapping): Namelist blocks returned by ``f90nml.read``.
        allowed_blocks (collection of str): Every block name accepted by the
        application. Names do not include the leading ``&``.

    Returns:
        None.

    Raises:
        ConfigError: If the namelist contains one or more blocks that are not
        in ``allowed_blocks``.
    """
    unknown_blocks = []

    for block_name in blocks:
        if block_name not in allowed_blocks:
            unknown_blocks.append(block_name)

    if unknown_blocks:
        formatted_names = []
        for block_name in unknown_blocks:
            formatted_names.append(f"&{block_name}")

        unknown_text = ", ".join(formatted_names)
        raise ConfigError(f"Unknown namelist block(s): {unknown_text}")


def require_fields(block, required_fields, block_label):
    """Check that required fields are present and do not contain ``None``.

    Args:
        block (mapping): One namelist block returned by ``f90nml.read``.
        required_fields (list or tuple of str): Field names that must be
        present in ``block``.
        block_label (str): Human-readable block label used in error messages,
        for example ``"&compare"``.

    Returns:
        None.

    Raises:
        ConfigError: If one or more required fields are missing or contain
        ``None``.
    """
    missing_fields = []

    for field_name in required_fields:
        field_is_missing = field_name not in block
        field_is_none = not field_is_missing and block[field_name] is None

        if field_is_missing or field_is_none:
            missing_fields.append(field_name)

    if missing_fields:
        missing_text = ", ".join(missing_fields)
        raise ConfigError(
            f"Missing required field(s) in {block_label}: {missing_text}"
        )


def reject_unknown_fields(block, allowed_fields, block_label):
    """Reject fields that are not part of an application's block schema.

    Args:
        block (mapping): One namelist block returned by ``f90nml.read``.
        allowed_fields (list or tuple of str): Every field name accepted in
        ``block``.
        block_label (str): Human-readable block label used in error messages,
        for example ``"&compare"``.

    Returns:
        None.

    Raises:
        ConfigError: If ``block`` contains one or more fields that are not in
        ``allowed_fields``.
    """
    unknown_fields = []

    for field_name in block:
        if field_name not in allowed_fields:
            unknown_fields.append(field_name)

    if unknown_fields:
        unknown_text = ", ".join(unknown_fields)
        raise ConfigError(f"Unknown field(s) in {block_label}: {unknown_text}")


def validate_schema(blocks, schema):
    """Validate namelist block and field structure against an application schema.

    Each schema entry is keyed by a namelist block name. Its value is a
    mapping containing ``required``, ``required_fields``, and
    ``optional_fields``.

    Args:
        blocks (mapping): Namelist blocks returned by ``f90nml.read``.
        schema (mapping): Application-owned block and field schema. Required
        and optional field collections contain strings.

    Returns:
        None.

    Raises:
        ConfigError: If required blocks or fields are missing, or if the
        namelist contains unknown blocks or fields.
    """

    # Step 1: collect and check the blocks marked as required by the schema.
    required_blocks = []
    for block_name, block_schema in schema.items():
        if block_schema["required"]:
            required_blocks.append(block_name)

    require_blocks(
        blocks=blocks,
        required_blocks=required_blocks,
    )

    # Step 2: reject block names that are not defined by the application.
    allowed_blocks = schema.keys()
    reject_unknown_blocks(
        blocks=blocks,
        allowed_blocks=allowed_blocks,
    )

    # Step 3: validate the fields in every block that was supplied. Optional
    # blocks that are absent require no further processing.
    for block_name, block_schema in schema.items():
        if block_name not in blocks:
            continue

        block = blocks[block_name]
        required_fields = block_schema["required_fields"]
        optional_fields = block_schema["optional_fields"]
        block_label = f"&{block_name}"

        allowed_fields = []
        for field_name in required_fields:
            allowed_fields.append(field_name)
        for field_name in optional_fields:
            allowed_fields.append(field_name)

        require_fields(
            block=block,
            required_fields=required_fields,
            block_label=block_label,
        )
        reject_unknown_fields(
            block=block,
            allowed_fields=allowed_fields,
            block_label=block_label,
        )


def normalize_string(value):
    """Normalize a selector string without hiding an invalid value type.

    Args:
        value (object): Value to normalize. String values are stripped of
        surrounding whitespace and converted to lowercase.

    Returns:
        object: The normalized string, or the original value when it is not a
        string. Preserving invalid types allows a later validation function to
        report them.
    """
    if isinstance(value, str):
        return value.strip().lower()
    return value


def normalize_path(value, config_path, field_label):
    """Resolve a configured path relative to its namelist file.

    Args:
        value (str or Path): Path value read from the configuration.
        config_path (str or Path): Path to the namelist file containing the
        configured value.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        Path: An absolute, normalized path. Relative values are resolved from
        the directory containing ``config_path``.

    Raises:
        ConfigError: If ``value`` is not a string or ``Path``, or represents an
        empty path.
    """
    if not isinstance(value, (str, Path)):
        raise ConfigError(f"{field_label} must be a string or Path.")

    path_text = str(value).strip()
    if not path_text:
        raise ConfigError(f"{field_label} must be a nonempty path.")

    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path

    return path.resolve()


def as_list(value):
    """Represent ``None``, a scalar, or a sequence consistently as a list.

    Args:
        value (object): Configuration value to normalize. It may be ``None``,
        a scalar, a list, or a tuple.

    Returns:
        list: An empty list for ``None``, a copy of a list, a tuple converted
        to a list, or a one-element list containing a scalar value.
    """
    if value is None:
        return []
    if isinstance(value, list):
        return value.copy()
    if isinstance(value, tuple):
        return list(value)
    return [value]


def require_string(value, field_label):
    """Require a nonempty string.

    Args:
        value (object): Configuration value to validate.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        str: The validated string with surrounding whitespace removed. Letter
        case is preserved.

    Raises:
        ConfigError: If ``value`` is not a string or is empty after stripping
        surrounding whitespace.
    """
    if not isinstance(value, str):
        raise ConfigError(f"{field_label} must be a string.")

    stripped_value = value.strip()
    if not stripped_value:
        raise ConfigError(f"{field_label} must be a nonempty string.")

    return stripped_value


def require_boolean(value, field_label):
    """Require a logical configuration value.

    Args:
        value (object): Configuration value to validate.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        bool: The validated logical value.

    Raises:
        ConfigError: If ``value`` is not a Python ``bool``.
    """
    if not isinstance(value, bool):
        raise ConfigError(f"{field_label} must be .true. or .false.")
    return value


def require_integer(value, field_label):
    """Require an integer configuration value and reject logical values.

    Args:
        value (object): Configuration value to validate.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        int: The validated integer value.

    Raises:
        ConfigError: If ``value`` is not an integer or is a Python ``bool``.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"{field_label} must be an integer.")
    return value


def require_real(value, field_label):
    """Require a real configuration value.

    Args:
        value (object): Configuration value to validate. Integers and
        floating-point values are accepted, but logical values are not.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        float: The validated value converted to a Python ``float``.

    Raises:
        ConfigError: If ``value`` is not an integer or floating-point value,
        or is a Python ``bool``.
    """
    valid_number = isinstance(value, (int, float)) and not isinstance(value, bool)
    if not valid_number:
        raise ConfigError(f"{field_label} must be a real number.")
    return float(value)


def require_choice(value, supported_values, field_label):
    """Require a value to belong to an application-supplied collection.

    Args:
        value (object): Configuration value to validate.
        supported_values (list, tuple, or set): Values accepted by the
        application.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        object: The validated value, unchanged.

    Raises:
        ConfigError: If ``value`` is not present in ``supported_values``.
    """
    if value not in supported_values:
        formatted_values = []
        for supported_value in supported_values:
            formatted_values.append(str(supported_value))

        supported_text = ", ".join(formatted_values)
        raise ConfigError(
            f"Unsupported {field_label} '{value}'. "
            f"Supported values: {supported_text}"
        )

    return value


def require_existing_file(path, field_label):
    """Require a path to identify an existing regular file.

    Args:
        path (str or Path): File path to validate.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        Path: The validated path represented as a ``Path`` object.

    Raises:
        ConfigError: If ``path`` does not identify an existing regular file.
    """
    path = Path(path)
    if not path.is_file():
        raise ConfigError(
            f"{field_label} does not identify an existing file: {path}"
        )
    return path


def require_existing_directory(path, field_label):
    """Require a path to identify an existing directory.

    Args:
        path (str or Path): Directory path to validate.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        Path: The validated path represented as a ``Path`` object.

    Raises:
        ConfigError: If ``path`` does not identify an existing directory.
    """
    path = Path(path)
    if not path.is_dir():
        raise ConfigError(
            f"{field_label} does not identify an existing directory: {path}"
        )
    return path
