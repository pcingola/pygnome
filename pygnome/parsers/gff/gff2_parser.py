"""
Parser for GFF2 format files.
"""

import re

from .base_parser import BaseParser
from .record import GffRecord
from pygnome.genomics.strand import Strand


class Gff2Parser(BaseParser):
    """Parser for GFF2 format files."""
    
    def _parse_line(self, line: str) -> GffRecord | None:
        """Parse a single line from a GFF2 file."""
        fields = line.strip().split('\t')
        
        # GFF2 has 9 fields
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
        """Parse GFF2 attributes into a dictionary."""
        attributes = {}
        
        # Handle empty attribute string
        if attr_string == '.' or not attr_string:
            return attributes
        
        # GFF2 attributes format: Class name [; attribute value]*
        # Example: "Transcript B0273.1; Note "Zn-Finger"; Alias MFX"
        
        # First, extract the class and name
        match = re.match(r'^(\S+)\s+(.+?)(?:;|$)', attr_string)
        if match:
            class_name, value = match.groups()
            attributes['class'] = class_name
            attributes['name'] = value.strip('"')
            
            # Remove the class/name part from the string
            attr_string = attr_string[match.end():]
        
        # Parse remaining attributes
        # Look for patterns like: Note "value" or Alias value
        for match in re.finditer(r';\s*(\w+)\s+(?:"([^"]+)"|([^;]+))(?:;|$)', attr_string):
            key, quoted_value, unquoted_value = match.groups()
            value = quoted_value if quoted_value is not None else unquoted_value
            
            # Handle multiple values for the same key (e.g., multiple Notes)
            if key in attributes:
                if isinstance(attributes[key], list):
                    attributes[key].append(value)
                else:
                    attributes[key] = [attributes[key], value]
            else:
                attributes[key] = value
        
        return attributes