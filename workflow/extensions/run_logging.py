"""Timestamped progress messages; preserve raw scientific program output formats."""
from datetime import datetime, timezone


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def log(message, *, file=None):
    print(f'[{timestamp()}] {message}', file=file, flush=True)
