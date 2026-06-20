# Command-line interface

APEX ships a headless CLI (`apex`, or `python -m apex`) that runs without a
display — for servers, CI, batch processing, and reproducible analysis. It never
imports PyQt5 at startup, so it works on a headless core install.

## Commands

### `apex doctor`

Diagnoses the runtime: Python version, required/optional dependencies, external
WCS solvers (ASTAP, astrometry.net), and the config file. Add `--network` to
also probe Gaia/SIMBAD reachability.

```bash
apex doctor
apex doctor --network
```

Exit code is non-zero only on a **blocking** problem (e.g. a missing required
dependency); optional gaps (no GUI, no external solver) are warnings.

### `apex config`

```bash
apex config init      # create parameters.toml from the bundled example
apex config path      # print the resolved parameters.toml location
apex config show      # print the active parameters.toml
```

### `apex run`

Runs the shared headless pipeline (Steps 1–7).

```bash
apex run --mode cmd --dry-run            # preview the plan
apex run --mode cmd --steps 1-7          # run all shared steps
apex run --mode cmd --steps 4-7 --force  # re-run a subset, ignoring cached output
apex run --mode lc  --result-dir D:\out  # override the output directory
```

Completed-step detection means a step whose outputs already exist is skipped
(use `--force` to override). Each run writes a `pipeline_run.json` manifest into
the result directory.

### `apex export`

Converts APEX light-curve output to community submission formats.

```bash
apex export --format aavso    --input lc.csv --output report.txt --obscode ABC --target "NGC6811-V1"
apex export --format exoclock --input lc.csv --output exoclock.txt
apex export --format exofop   --input lc.csv --output tfop.csv
```

### `apex gui`

Launches the desktop app (requires the `[gui]` extra).

```bash
apex gui            # launcher
apex gui --mode cmd # open a mode directly
```

### `apex version`

Prints the APEX and Python versions.
