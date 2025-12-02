from rich.console import Console
from rich.table import Table
from rich import box
import pyfiglet
import sys
import io
import logging
from contextlib import contextmanager, redirect_stderr, redirect_stdout

# Exported, shared console for the app
console = Console()


def print_separator(title: str | None = None) -> None:
    """Print a simple, consistent separator using hyphens."""
    line = "-" * 100
    if title and title.strip():
        console.print(line, style="dim")
        console.print(title.strip(), style="bold white")
        console.print(line, style="dim")
    else:
        console.print(line, style="dim")


def print_header() -> None:
    """Render the application header using pyfiglet and Rich (single color, readable)."""
    title = pyfiglet.figlet_format("edna2obis", font="standard")
    console.print(title, style="navy_blue")
    console.print("eDNA to OBIS/GBIF Publisher (v3) - NOAA Omics", style="bold navy_blue")
    print_separator()


def print_usage() -> None:
    """Show primary usage modes in a clean, scannable table."""
    table = Table(title="Usage", box=box.SIMPLE, show_lines=False, header_style="bold")
    table.add_column("Mode", style="navy_blue", no_wrap=True)
    table.add_column("Command", style="white")
    table.add_column("Description", style="dim")

    # Mode 1: Main pipeline (preferred)
    table.add_row(
        "Mode 1",
        "python main.py",
        "Main pipeline: Occurrence Core, DNA Derived Extension, optional eMoF, optional EML, "
        "Taxa Assignment INFO, and HTML Report",
    )

    # Mode 2: Taxonomy-only via CLI wrapper (do not call the script directly)
    table.add_row(
        "Mode 2",
        "python taxassign.py --help",
        "Taxonomy-only: Assign taxonomy via WoRMS/GBIF to a column of names",
    )

    console.print(table)
    console.print(
        "For flags and examples, run `python taxassign.py --help`. See README for advanced configuration.",
        style="dim",
    )


@contextmanager
def silence_output(min_log_level: int = logging.WARNING):
    """
    Suppress stderr (logging, warnings) and temporarily raise logging level.
    Leaves stdout untouched so spinners and normal console output remain visible.
    """
    fake_err = io.StringIO()
    root_logger = logging.getLogger()
    old_level = root_logger.level
    try:
        root_logger.setLevel(min_log_level)
        with redirect_stderr(fake_err):
            yield
    finally:
        root_logger.setLevel(old_level)


@contextmanager
def silence_stdout():
    """Temporarily redirect stdout to keep Rich status line from being pushed up."""
    fake_out = io.StringIO()
    with redirect_stdout(fake_out):
        yield


@contextmanager
def silence_stdouterr(min_log_level: int = logging.WARNING):
    """
    Suppress both stdout and stderr and raise logging level.
    Use when third-party code prints to stdout and would disrupt the spinner.
    """
    fake_out = io.StringIO()
    fake_err = io.StringIO()
    root_logger = logging.getLogger()
    old_level = root_logger.level
    try:
        root_logger.setLevel(min_log_level)
        with redirect_stdout(fake_out), redirect_stderr(fake_err):
            yield
    finally:
        root_logger.setLevel(old_level)

