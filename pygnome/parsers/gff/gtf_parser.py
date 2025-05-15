"""
Parser for GTF format files.
"""

import re

from .base_parser import BaseParser
from .record import GffRecord
from pygnome.genomics.strand import Strand


class GtfParser(BaseParser):
    """Parser for GTF format files."""
    
    def _parse_line(self, line: str) -> GffRecord | None:
        """Parse a single line from a GTF file."""
        fields = line.strip().split('\t')
        
        # GTF has 9 fields
        if len(fields) != 9:
            return None
        
        # Extract basic fields
        seqid, source, type_, start, end, score, strand, phase, attr_string = fields
        
        # Parse attributes
        attributes = self._parse_attributes(attr_string)
        
        # Create and return record
        return GffRecord(
            id=attributes.get("ID", f"{type_}_{seqid}_{start}_{end}"),
            chrom=seqid,
            source=source,
            type=type_,
            start=start,
            end=end,
            score=score,
            strand=strand,
            phase=phase,
            attributes=attributes
        )
    
    def _parse_attributes(self, attr_string: str) -> dict:
        """Parse GTF attributes into a dictionary."""
        attributes = {}
        
        # Handle empty attribute string
        if attr_string == '.' or not attr_string:
            return attributes
        
        # GTF attributes format: tag "value"; tag "value";
        # Extract all key-value pairs
        for match in re.finditer(r'(\S+)\s+"([^"]+)"\s*;', attr_string):
            key, value = match.groups()
            attributes[key] = value
        
        return attributes