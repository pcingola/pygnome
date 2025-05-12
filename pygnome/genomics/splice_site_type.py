"""Splice site type enumeration."""

from enum import Enum


class SpliceSiteType(str, Enum):
    """Type of splice site."""
    DONOR = "donor"      # 5' splice site (GT)
    ACCEPTOR = "acceptor"  # 3' splice site (AG)