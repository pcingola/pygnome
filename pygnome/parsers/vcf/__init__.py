"""
VCF Parser module for pygnome.

This module provides classes for parsing and working with VCF (Variant Call Format) files.
"""

from pygnome.parsers.vcf.vcf_reader import VcfReader
from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.variant_factory import VariantFactory

__all__ = [
    "VcfReader",
    "VcfHeader",
    "VcfRecord",
    "VariantFactory",
]