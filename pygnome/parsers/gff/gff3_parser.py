"""
Parser for GFF3 format files.
"""

from urllib.parse import unquote

from .base_parser import BaseParser
from .record import GffRecord
from pygnome.genomics.strand import Strand


class Gff3Parser(BaseParser):
    """Parser for GFF3 format files."""
    
    def _parse_line(self, line: str) -> GffRecord | None:
        """Parse a single line from a GFF3 file."""
        fields = line.strip().split('\t')
        
        # GFF3 has 9 fields
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
        """Parse GFF3 attributes into a dictionary."""
        attributes = {}
        
        # Handle empty attribute string
        if attr_string == '.' or not attr_string:
            return attributes
        
        # GFF3 attributes format: tag=value;tag=value
        for attr in attr_string.split(';'):
            if not attr:
                continue
                
            # Split on first equals sign
            parts = attr.strip().split('=', 1)
            if len(parts) != 2:
                continue
                
            key, value = parts
            
            # URL-decode the value
            value = unquote(value)
            
            # Handle attributes that can have multiple values (comma-separated)
            if key in ('Parent', 'Alias', 'Note', 'Dbxref', 'Ontology_term'):
                values = value.split(',')
                if key in attributes:
                    if isinstance(attributes[key], list):
                        attributes[key].extend(values)
                    else:
                        attributes[key] = [attributes[key]] + values
                else:
                    attributes[key] = values if len(values) > 1 else values[0]
            else:
                attributes[key] = value
        
        return attributes