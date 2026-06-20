# APEX Documentation

## User and Operator Guides

- [Project README](../README.md): installation, workflow, WCS overview
- [Configuration Guide](configuration.md): `parameters.toml` ownership and
  major setting groups
- [Windows Deployment](../deploy/README.md): installer and portable release
- [Benchmark Guide](../benchmark/README.md): artificial-star and IRAF checks
- [Validation Guide](../validation/README.md): synthetic and real-data checks

## Architecture and Design

- [Architecture](../ARCHITECTURE.md): modules, workflow dispatch, data flow
- [Cache Manager Design](cache-manager-design.md): cache ownership and
  invalidation policy
- [Cache Inventory](cache-inventory.md): cache-like paths in the current code
- [PSF Photometry Analysis](psf-photometry-analysis.md): algorithm and observed
  limitations

## Generated or Planning Documents

- [Parameter Inventory](parameter-inventory.md): generated TOML/runtime mapping
  inventory
- [Config/Cache/UI Refactor Plan](config-cache-ui-refactor-plan.md): migration
  plan, not a current behavior specification

When documents disagree, use these sources of truth:

1. `apex/gui/main_window.py` for step dispatch
2. `apex/utils/step_paths*.py` for output paths
3. `parameters.example.toml` and mode parameter maps for defaults
4. `deploy/build_release.bat` for release behavior
