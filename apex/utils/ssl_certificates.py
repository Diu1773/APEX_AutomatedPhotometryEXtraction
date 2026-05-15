"""HTTPS certificate setup for packaged astronomy web queries."""

from __future__ import annotations

import os
from pathlib import Path


def configure_ssl_certificates() -> tuple[bool, str]:
    """Point Python HTTPS clients at certifi's CA bundle when available."""
    try:
        import certifi
    except Exception as exc:
        return False, f"certifi is not importable: {type(exc).__name__}: {exc}"

    try:
        ca_path = Path(certifi.where())
    except Exception as exc:
        return False, f"certifi.where() failed: {type(exc).__name__}: {exc}"

    if not ca_path.is_file():
        return False, f"certifi CA bundle not found: {ca_path}"

    for key in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE"):
        current = os.environ.get(key)
        if not current or not Path(current).is_file():
            os.environ[key] = str(ca_path)
    return True, f"certifi CA bundle: {ca_path}"


def ssl_error_reason(value: object) -> str | None:
    text = str(value or "").lower()
    if (
        "certificate_verify_failed" in text
        or "certificate verify failed" in text
        or "unable to get local issuer certificate" in text
        or "sslcertverificationerror" in text
    ):
        return "ssl_cert_verify_failed"
    return None
