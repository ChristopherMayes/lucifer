"""The tables the package's objects show in Jupyter, through their _repr_html_."""
from __future__ import annotations

import html


def table(title, rows):
    """An HTML table with a caption, ``title``, and one row per (name, value) pair."""
    cells = "".join(
        f"<tr><th style='text-align:left'>{html.escape(str(k))}</th>"
        f"<td style='text-align:left'>{html.escape(str(v))}</td></tr>" for k, v in rows)
    return (f"<table><caption style='text-align:left'><b>{html.escape(title)}</b></caption>"
            f"{cells}</table>")


def short(value, limit=60):
    """A value as text, cut at ``limit`` characters."""
    text = repr(value) if isinstance(value, str) else str(value)
    return text if len(text) <= limit else text[:limit - 3] + "..."
