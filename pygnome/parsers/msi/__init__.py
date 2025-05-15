"""
MSI-sites parser module.

This module provides classes for parsing MSI (Microsatellite Instability) sites files,
which contain information about microsatellite locations in the genome.
"""

from .msi_site_record import MsiSiteRecord
from .msi_sites_reader import MsiSitesReader

__all__ = ["MsiSiteRecord", "MsiSitesReader"]