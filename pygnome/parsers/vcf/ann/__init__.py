"""
VCF ANN field parser module.

This module provides classes for parsing the ANN field in VCF files,
which contains variant annotation information.
"""

from .ann_parser import AnnParser
from .variant_annotation import VariantAnnotation
from .enums import (
    AnnotationImpact,
    FeatureType,
    ErrorWarningType
)

__all__ = [
    "AnnParser",
    "VariantAnnotation",
    "AnnotationImpact",
    "FeatureType",
    "ErrorWarningType"
]