"""MkDocs hooks for the tissue_simulator docs site.

Silences griffe "No type or annotation for parameter ..." warnings emitted
when mkdocstrings introspects a handful of older callbacks in the
tissue_simulator source tree. The docs-site scaffolding is not allowed to
modify the tissue_simulator source, so we filter these specific log records
at the griffe logger level instead. Real regressions (broken cross-refs,
missing pages, etc.) still surface through mkdocs --strict.
"""

from __future__ import annotations

import logging
import re

_GRIFFE_NOISE = re.compile(
    r"No type or annotation for (?:parameter|returned value) "
)


class _GriffeAnnotationFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:  # noqa: D401
        msg = record.getMessage()
        if _GRIFFE_NOISE.search(msg):
            return False
        return True


def on_startup(command, dirty):  # noqa: ARG001 - mkdocs hook signature
    """Attach the griffe-annotation filter as early as possible.

    Empirically the warnings come from the ``mkdocs.plugins.griffe``
    logger; we also attach to a couple of nearby names so this stays
    robust if mkdocstrings ever renames its wrapper logger.
    """
    flt = _GriffeAnnotationFilter()
    for name in (
        "mkdocs.plugins.griffe",
        "mkdocs.plugins.mkdocstrings",
        "griffe",
        "mkdocs",
    ):
        logging.getLogger(name).addFilter(flt)
