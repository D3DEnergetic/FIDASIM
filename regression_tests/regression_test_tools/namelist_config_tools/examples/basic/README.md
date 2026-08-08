# Basic configuration example

This directory demonstrates how an application-specific `config.py` can use
the generic `namelist_config_tools` package directly.

From the FIDASIM repository root, run:

```bash
python3 \
  regression_tests/regression_test_tools/namelist_config_tools/examples/basic/show_config.py \
  regression_tests/regression_test_tools/namelist_config_tools/examples/basic/input_config.nml
```

The example validates the namelist and prints the canonical configuration
keys followed by their normalized values. Relative paths work independently
of the directory from which the command is launched.

The files have separate responsibilities:

| File | Responsibility |
| --- | --- |
| `input_config.nml` | User-facing Fortran namelist. |
| `config.py` | Application schema, semantic checks, defaults, and canonical return dictionary. |
| `show_config.py` | Thin command-line consumer and `ConfigError` boundary. |
| `input_data.txt` | Existing input used to demonstrate path validation. |

Try introducing an unknown field, choosing an unsupported mode, making
`n_steps` negative, or removing `&output` while `write_output = .true.` to see
the resulting configuration errors.
