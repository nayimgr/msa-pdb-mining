"""Output renderers: machine-readable data, static image, interactive HTML."""

from __future__ import annotations

from .data import write_data_outputs
from .html import write_html
from .image import write_image

__all__ = ["write_data_outputs", "write_html", "write_image"]
