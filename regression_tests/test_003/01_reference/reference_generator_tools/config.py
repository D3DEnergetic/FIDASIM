"""Read and validate the Test 003 reference-validation configuration."""

from regression_test_tools import (
    normalize_path,
    read_namelist,
    require_existing_file,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "reference": {
        "required": True,
        "required_fields": ["input_distribution_config"],
        "optional_fields": ["comment"],
    },
}


def read_config(config_filename):
    """Return the resolved Test 002 Stage 2 configuration to validate."""
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)
    reference_block = blocks["reference"]

    distribution_config = normalize_path(
        value=require_string(
            value=reference_block["input_distribution_config"],
            field_label="input_distribution_config",
        ),
        config_path=config_path,
        field_label="input_distribution_config",
    )
    require_existing_file(
        path=distribution_config,
        field_label="input_distribution_config",
    )

    comment = ""
    if "comment" in reference_block:
        comment = require_string(
            value=reference_block["comment"],
            field_label="comment",
        )

    return {
        "reference": {
            "comment": comment,
            "input_distribution_config": distribution_config,
        }
    }
