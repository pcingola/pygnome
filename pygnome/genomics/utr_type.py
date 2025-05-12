"""UTR type enumeration."""

from enum import Enum


class UTRType(str, Enum):
    """Type of UTR (untranslated region)."""
    FIVE_PRIME = "5'"
    THREE_PRIME = "3'"