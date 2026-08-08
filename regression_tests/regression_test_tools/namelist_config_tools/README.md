# Generic namelist configuration tools

`namelist_config_tools` provides common, application-independent operations
for configuration workflows based on Fortran namelists. It reads namelists
with `f90nml`, validates their declared structure, normalizes common value
forms, and reports configuration mistakes through one exception type.

The module is intentionally not a complete configuration system. A local
`config.py` remains responsible for defining the application's schema,
physical constraints, conditional requirements, defaults, and canonical
return structure.

The reusable pattern is therefore:

```text
input_config.nml
       |
       v
namelist_config_tools
  - read the namelist
  - validate blocks and field names
  - normalize common types and paths
       |
       v
application config.py
  - enforce ranges and coupled conditions
  - add defaults
  - return one canonical dictionary
       |
       v
calculation or plotting workflow
```

## Dependency and import

The only third-party dependency of `config.py` is
[`f90nml`](https://pypi.org/project/f90nml/):

```bash
python3 -m pip install f90nml
```

Within FIDASIM, the outer `regression_test_tools` package is a compatibility
facade. Regression workflows use only this established public interface:

```python
from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    require_integer,
    require_string,
    validate_schema,
)
```

FIDASIM code does not import the nested implementation directly. In another
project, copy the complete `namelist_config_tools` directory into an
importable source directory and use its direct interface:

```python
from namelist_config_tools import (
    ConfigError,
    read_namelist,
    validate_schema,
)
```

The copied directory does not yet provide packaging metadata for installation
with `pip`; place its parent directory on `PYTHONPATH` or within the other
project's importable source tree.

## The two validation layers

### 1. Generic structural validation

The shared module answers questions that are independent of the application:

- Does the namelist file exist and parse successfully?
- Are all required namelist blocks present?
- Are any unknown blocks or fields present?
- Are required fields present and not unset (`None`) at the namelist level?
- Is a value a string, logical, integer, or real number?
- Is a selector one of the application-supplied choices?
- Is a referenced file or directory present?
- How should a relative path be resolved from the configuration file?

### 2. Application-specific semantic validation

The local `config.py` answers questions that require domain knowledge:

- Must an integer be positive or fit within a Fortran integer kind?
- Must a real value be finite and lie within a physical range?
- Is one block required only when a selector or logical is active?
- Must two arrays have equal lengths?
- Is a filename required to have a particular extension?
- Which defaults should be inserted into the canonical configuration?
- Which derived paths or values should downstream code receive?

Keeping these layers separate makes the shared module reusable without
embedding Test 002, Test 003, Test 004, FIDASIM, or physics-specific policy in
it.

## Public API

All supported public names are exported by `namelist_config_tools/__init__.py`.
FIDASIM's outer `regression_test_tools/__init__.py` re-exports the same API.

| Function or type | Purpose |
| --- | --- |
| `ConfigError` | Specialized `ValueError` used for expected configuration failures. Command-line entry points can catch it without hiding unrelated programming errors. |
| `read_namelist(config_path)` | Resolve and read a namelist. Returns `(resolved_config_path, blocks)`. |
| `validate_schema(blocks, schema)` | Validate required blocks, required fields, and reject unknown blocks and fields. |
| `require_blocks(blocks, required_blocks)` | Require selected blocks when the condition is determined dynamically. |
| `require_fields(block, required_fields, block_label)` | Require selected fields when a complete static schema is unnecessary. |
| `reject_unknown_blocks(blocks, allowed_blocks)` | Reject block names outside an application-owned collection. |
| `reject_unknown_fields(block, allowed_fields, block_label)` | Reject field names outside an application-owned collection. |
| `normalize_string(value)` | Strip and lowercase strings while leaving invalid non-string values unchanged for later validation. |
| `normalize_path(value, config_path, field_label)` | Resolve a path relative to the namelist that contains it and return an absolute `Path`. It does not require the target to exist. |
| `as_list(value)` | Normalize `None`, scalars, lists, and tuples to a new list. |
| `require_string(value, field_label)` | Require a nonempty string and return it stripped, with letter case preserved. |
| `require_boolean(value, field_label)` | Require a namelist logical value and return a Python `bool`. |
| `require_integer(value, field_label)` | Require a Python integer and reject logical values. |
| `require_real(value, field_label)` | Accept an integer or float, reject logical values, and return a Python `float`. |
| `require_choice(value, supported_values, field_label)` | Require membership in an application-supplied collection. |
| `require_existing_file(path, field_label)` | Require an existing regular file and return it as a `Path`. |
| `require_existing_directory(path, field_label)` | Require an existing directory and return it as a `Path`. |
| `print_config(config)` | Print the keys of a nested canonical configuration as a tree. Values are deliberately omitted. |

## Schema format

`validate_schema` expects a dictionary keyed by namelist block name:

```python
CONFIG_SCHEMA = {
    "run": {
        "required": True,
        "required_fields": ["input_filename", "n_steps"],
        "optional_fields": ["comment"],
    },
    "output": {
        "required": False,
        "required_fields": ["output_directory"],
        "optional_fields": [],
    },
}
```

This schema means:

- `&run` must be present;
- `input_filename` and `n_steps` must occur in `&run`;
- `comment` may occur in `&run`;
- `&output` may be absent;
- when `&output` is present, `output_directory` is required; and
- every unlisted block or field is rejected as a likely typo.

The schema checks structure only. It does not validate value types, ranges,
units, file extensions, or relationships between settings.

## Recommended local `config.py` workflow

A local reader normally performs these steps in order:

```python
def read_config(config_filename):
    # 1. Read the namelist and retain its resolved path.
    config_path, blocks = read_namelist(config_path=config_filename)

    # 2. Reject missing and unknown interface elements.
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    # 3. Select the validated blocks.
    run_block = blocks["run"]

    # 4. Validate and normalize each value explicitly.
    mode = require_choice(
        value=normalize_string(run_block["mode"]),
        supported_values=["fast", "accurate"],
        field_label="mode",
    )
    n_steps = require_integer(run_block["n_steps"], "n_steps")

    # 5. Add semantic checks owned by this application.
    if n_steps < 1:
        raise ConfigError("n_steps must be positive.")

    # 6. Resolve paths from the file in which they were written.
    input_filename = normalize_path(
        value=require_string(
            run_block["input_filename"],
            "input_filename",
        ),
        config_path=config_path,
        field_label="input_filename",
    )
    require_existing_file(input_filename, "input_filename")

    # 7. Return a stable canonical API for downstream code.
    return {
        "run": {
            "mode": mode,
            "n_steps": n_steps,
            "input_filename": input_filename,
        }
    }
```

Downstream calculations should consume only this canonical dictionary. They
should not repeatedly interpret raw namelist values.

## Important behavior

### Relative paths belong to their namelist

`normalize_path` resolves a relative value from the directory containing the
configuration file, not from the caller's current working directory. This
makes a configuration behave consistently when invoked from different
locations.

```python
input_filename = normalize_path(
    value="data/input.h5",
    config_path="/project/cases/case_A.nml",
    field_label="input_filename",
)

# /project/cases/data/input.h5
```

Use `require_existing_file` or `require_existing_directory` after
normalization only when the path is an input that must already exist. Output
paths normally should not be required to exist.

### Normalization does not replace validation

`normalize_string` intentionally returns a non-string unchanged. This allows
`require_string`, `require_choice`, or a local validator to produce the
appropriate error instead of silently coercing a bad value.

Likewise, `as_list` standardizes cardinality but does not validate the list
elements. Validate every normalized element explicitly.

### Logical values are not integers

Python treats `bool` as a subclass of `int`. `require_integer` and
`require_real` explicitly reject logical values so `.true.` cannot be
accidentally accepted as `1`.

### Required is different from active

A schema describes what is structurally allowed. Conditional activation is a
semantic rule. For example:

```python
if write_output:
    require_blocks(blocks=blocks, required_blocks=["output"])
```

Test 004 uses this pattern when plotting or saving activates blocks that are
otherwise optional.

### Keep error labels user-facing

The `field_label` argument appears in error messages. Prefer explicit calls:

```python
seed = require_integer(
    value=run_block["seed"],
    field_label="seed",
)
```

The label describes the configuration concept; it is not necessarily a Python
variable name or filename.

## Error handling at a command-line boundary

Local readers raise `ConfigError`. A thin command-line entry point should
translate that expected failure into a concise message and nonzero exit code:

```python
try:
    config = read_config(arguments.config_path)
except ConfigError as error:
    raise SystemExit(f"Configuration error: {error}") from None
```

Avoid catching every `Exception`: unexpected programming errors should retain
their traceback during development.

## Patterns used by Tests 002–004

The regression tests demonstrate progressively richer uses of the same
generic layer:

- **[Test 002 Stage 1](../../test_002/01_reference/reference_tools/config.py)**
  uses `as_list` to normalize scalar or repeated spatial locations, validates
  each element, normalizes selectors, and resolves the CQL3D input path.
- **[Test 002 Stage 2](../../test_002/02_run_test/conversion_tools/config.py)**
  uses a complete static schema, validates a required upstream configuration,
  applies numerical grid constraints, and returns a canonical
  block-structured dictionary.
- **[Test 003 Stage 2](../../test_003/02_run_test/normalize_config.py)** uses
  `require_blocks` and `require_fields` when it needs only a small contract
  from an upstream configuration owned by another stage.
- **[Test 003 Stage 3](../../test_003/03_compare/comparison_tools/config.py)**
  combines a static comparison schema with file discovery and validation of
  generated input collections.
- **[Test 004](../../test_004/test_004_tools/config.py)** defines shared and
  implementation-specific schemas. It makes inactive implementation blocks
  optional and uses dynamic `require_blocks` checks when plotting or output
  settings activate another block.
- **Test 004 Stage 2** converts the canonical nested dictionary into a flat,
  Fortran-readable namelist artifact. This keeps the user-facing interface in
  Python while allowing the numerical calculation to remain in Fortran.

These examples also show what should remain local: finite-value checks,
positive ranges, integer-kind limits, selector-dependent fields, safe output
rules, and derived filenames.

## Runnable example

[`examples/basic`](examples/basic/README.md) contains a complete minimal
application:

```text
examples/basic/
├── README.md
├── config.py
├── input_config.nml
├── input_data.txt
└── show_config.py
```

It demonstrates a static schema, scalar and list normalization, selector and
range validation, relative input and output paths, a conditionally required
output block, a canonical return dictionary, `print_config`, and CLI error
handling.

## Reusing the module in another project

For a small internal project:

1. Copy `namelist_config_tools/` into the project's importable source tree.
2. Install `f90nml` in the project's Python environment.
3. Keep the shared `config.py` unchanged.
4. Create one application-local `config.py` containing its schema and
   semantic rules.
5. Make calculations consume only the canonical dictionary returned by that
   local reader.
6. Add valid and invalid configuration cases to the project's tests.

If the tool is shared by several repositories, the next natural step is to
extract it into a small independently versioned Python package with packaging
metadata and unit tests. Its current API is already separated well enough for
that extraction, but this repository intentionally keeps it local for now.
