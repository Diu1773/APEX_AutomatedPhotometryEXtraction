"""Enable ``python -m apex`` as an alias for the headless CLI."""

from apex.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
