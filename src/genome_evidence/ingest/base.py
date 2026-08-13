"""Shared ingestion types."""

from enum import StrEnum


class ParseMode(StrEnum):
    STRICT = "strict"
    LENIENT = "lenient"
