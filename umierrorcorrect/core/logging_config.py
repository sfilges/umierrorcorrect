#!/usr/bin/env python3
"""Logging configuration using loguru for umierrorcorrect."""

import contextlib
import sys
from datetime import datetime
from pathlib import Path
from typing import Literal

from loguru import logger

# Track file handler ID so we can avoid duplicates
_file_handler_id: int | None = None

# Default log format
LOG_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
    "<level>{message}</level>"
)

# Simple format for console output
SIMPLE_FORMAT = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>"


def setup_logging(
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO",
    log_file: str | Path | None = None,
    simple_format: bool = True,
    rotation: str = "10 MB",
    retention: str = "7 days",
) -> None:
    """Configure loguru logging for the application.

    Args:
        level: Minimum log level to display.
        log_file: Optional path to log file. If None, logs only to console.
        simple_format: Use simple format for console (default True).
        rotation: When to rotate log file (e.g., "10 MB", "1 day").
        retention: How long to keep old log files.
    """
    # Remove default handler
    logger.remove()

    # Console handler with colored output
    log_format = SIMPLE_FORMAT if simple_format else LOG_FORMAT
    logger.add(
        sys.stderr,
        format=log_format,
        level=level,
        colorize=True,
    )

    # File handler if specified
    if log_file:
        log_path = Path(log_file)
        logger.add(
            str(log_path),
            format=LOG_FORMAT,
            level=level,
            rotation=rotation,
            retention=retention,
            compression="gz",
        )


def get_logger(name: str | None = None):
    """Get a logger instance.

    Args:
        name: Optional name for the logger context.

    Returns:
        Configured logger instance.
    """
    if name:
        return logger.bind(name=name)
    return logger


# Intercept standard logging to redirect to loguru
class InterceptHandler:
    """Handler to intercept standard logging and redirect to loguru."""

    def __init__(self, level: int = 0) -> None:
        self.level = level

    def emit(self, record) -> None:
        """Emit a log record."""
        # Get corresponding Loguru level if it exists
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find caller from where originated the logged message
        frame, depth = sys._getframe(6), 6
        while frame and frame.f_code.co_filename == __file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())


def redirect_standard_logging() -> None:
    """Redirect standard logging module to use loguru."""
    import logging

    logging.basicConfig(handlers=[InterceptHandler()], level=0, force=True)


def get_log_path(output_dir: Path | str) -> Path:
    """Generate timestamped log file path.

    Args:
        output_dir: Directory where log file will be created.

    Returns:
        Path to the log file with format: simsen-cli_YYYYMMDD_HHMMSS.log
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(output_dir) / f"simsen-cli_{timestamp}.log"


def add_file_handler(
    log_path: Path | str,
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "DEBUG",
) -> int:
    """Add a file handler to the logger.

    This function ensures only one file handler is active at a time to avoid
    duplicate logs when multiple subcommands are run.

    Args:
        log_path: Path to the log file.
        level: Minimum log level for file logging.

    Returns:
        Handler ID that can be used to remove the handler later.
    """
    global _file_handler_id

    # Remove existing file handler if any
    if _file_handler_id is not None:
        with contextlib.suppress(ValueError):
            logger.remove(_file_handler_id)

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    _file_handler_id = logger.add(
        str(log_path),
        format=LOG_FORMAT,
        level=level,
        rotation="10 MB",
        retention="7 days",
        compression="gz",
    )

    return _file_handler_id


def log_subprocess_stderr(stderr: str | bytes | None, tool_name: str) -> None:
    """Log captured stderr from an external tool.

    Args:
        stderr: Stderr output from subprocess (string or bytes).
        tool_name: Name of the external tool for log context.
    """
    if stderr is None:
        return

    # Convert bytes to string if needed
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf-8", errors="replace")

    stderr = stderr.strip()
    if stderr:
        # Log each line separately for readability
        for line in stderr.splitlines():
            if line.strip():
                logger.debug(f"[{tool_name}] {line}")
