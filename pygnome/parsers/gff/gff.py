"""
Main interface for reading and working with GFF/GTF files.
"""

from pathlib import Path
from typing import Iterator, List

from .base_parser import BaseParser
from .gff2_parser import Gff2Parser
from .gff3_parser import Gff3Parser
from .gtf_parser import GtfParser
from .feature_hierarchy import FeatureHierarchy
from .record import GffRecord
from .format import Format


class Gff:
    """
    Main interface for reading and working with GFF/GTF files.
    
    This class is implemented as a context manager for efficient file handling
    and provides methods for filtering and accessing features.
    """
    
    def __init__(self, file_path: str, format: Format | None = None):
        """
        Initialize the Gff context manager.
        
        Args:
            file_path: Path to the GFF/GTF file
            format: Optional format specification (Format.GFF2, Format.GFF3, Format.GTF)
                   If not provided, format will be auto-detected
        """
        self.file_path = file_path
        self.format = format
        self.parser = None
        self.file_handle = None
        self.records = []
        self._hierarchy = None
    
    def __enter__(self):
        """Enter the context manager."""
        # Auto-detect format if not specified
        if not self.format:
            self.format = self._detect_format()
        
        # Create appropriate parser
        if self.format == Format.GFF2:
            self.parser = Gff2Parser()
        elif self.format == Format.GFF3:
            self.parser = Gff3Parser()
        elif self.format == Format.GTF:
            self.parser = GtfParser()
        else:
            raise ValueError(f"Unsupported format: {self.format}")
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the context manager."""
        # Clean up resources
        self.file_handle = None
        
        # Don't suppress exceptions
        return False
    
    def __iter__(self) -> Iterator[GffRecord]:
        """Iterate over records in the file."""
        # Parse the file and yield records
        for record in self.parser.parse(self.file_path):
            self.records.append(record)
            yield record
    
    def _detect_format(self) -> Format:
        """Detect the format of the file."""
        return BaseParser.detect_format(self.file_path)
    
    def build_hierarchy(self) -> FeatureHierarchy:
        """
        Build and return a feature hierarchy for parent-child relationships.
        
        Returns:
            FeatureHierarchy object for navigating feature relationships
        """
        if self._hierarchy is None:
            self._hierarchy = FeatureHierarchy()
            self._hierarchy.build_from_records(self.records)
        
        return self._hierarchy
    
    def get_features_by_type(self, type_: str) -> List[GffRecord]:
        """
        Get all features of a specific type (e.g., 'gene', 'exon', 'CDS').
        
        Args:
            type_: Feature type to filter by
            
        Returns:
            List of matching Record objects
        """
        return [r for r in self.records if r.type == type_]
    
    def get_features_by_location(self, seqid: str, start: int, end: int) -> List[GffRecord]:
        """
        Get all features that overlap with the specified genomic region.
        
        Args:
            seqid: Reference sequence identifier
            start: 1-based start coordinate of the region
            end: End coordinate of the region (inclusive)
            
        Returns:
            List of features that overlap with the region
        """
        return [
            r for r in self.records
            if r.chrom == seqid and
               not (r.end < start or r.start > end)  # Overlap check
        ]
    
    def get_features_by_attribute(self, attr_name: str, attr_value: str | None = None) -> List[GffRecord]:
        """
        Get features based on their attributes.
        
        Args:
            attr_name: Name of the attribute to search for (e.g., 'ID', 'Name', 'gene_id')
            attr_value: If provided, only return features where the attribute equals this value
            
        Returns:
            List of features with matching attributes
        """
        if attr_value is None:
            # Just check for presence of attribute
            return [r for r in self.records if attr_name in r.attributes]
        else:
            # Check for specific value
            return [
                r for r in self.records 
                if attr_name in r.attributes and r.attributes[attr_name] == attr_value
            ]