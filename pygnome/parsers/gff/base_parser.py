"""
Base parser class for GFF/GTF parsers.
"""

from abc import ABC, abstractmethod
from pathlib import Path
import gzip
import re
from typing import Iterator

from .record import GffRecord
from ...genomics.strand import Strand
from .format import Format


class BaseParser(ABC):
    """
    Abstract base class for GFF/GTF parsers.
    
    This class defines the common interface for all format-specific parsers
    and provides shared functionality like file handling and format detection.
    """
    
    @abstractmethod
    def _parse_line(self, line: str) -> GffRecord | None:
        """Parse a single line from the file into a Record object."""
        pass
    
    @abstractmethod
    def _parse_attributes(self, attr_string: str) -> dict:
        """Parse the attributes field into a dictionary."""
        pass
    
    def _validate_record(self, record: GffRecord) -> bool:
        """Validate a parsed record."""
        # Basic validation: check required fields
        if not record.chrom or not record.type:
            return False
        
        # Check that start <= end
        if record.start > record.end:
            return False
            
        return True
    
    def parse(self, file_path: str) -> Iterator[GffRecord]:
        """
        Parse a GFF/GTF file and yield Record objects.
        
        Args:
            file_path: Path to the GFF/GTF file
            
        Yields:
            Record objects for each valid feature in the file
        """
        path = Path(file_path)
        
        # Check if file exists
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Open file (handle gzipped files)
        open_func = gzip.open if str(file_path).endswith(('.gz', '.gzip')) else open
        
        with open_func(file_path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                try:
                    record = self._parse_line(line)
                    if record and self._validate_record(record):
                        yield record
                except Exception as e:
                    # Log error but continue parsing
                    print(f"Error parsing line {line_num}: {line}")
                    print(f"Error: {str(e)}")
    
    @classmethod
    def detect_format(cls, file_path: str) -> Format:
        """
        Detect the format of a GFF/GTF file.
        
        Args:
            file_path: Path to the GFF/GTF file
            
        Returns:
            Format enum value: Format.GFF2, Format.GFF3, or Format.GTF
        """
        path = Path(file_path)
        
        # Check file extension first
        if path.suffix.lower() in ('.gtf', '.gtf.gz'):
            return Format.GTF
        elif path.suffix.lower() in ('.gff3', '.gff3.gz'):
            return Format.GFF3
        elif path.suffix.lower() in ('.gff', '.gff.gz'):
            # Need to check content to distinguish GFF2 vs GFF3
            pass
        
        # Open file (handle gzipped files)
        open_func = gzip.open if str(file_path).endswith(('.gz', '.gzip')) else open
        
        with open_func(file_path, 'rt') as f:
            # Check first few non-comment lines
            for _ in range(100):  # Check up to 100 lines
                line = f.readline()
                if not line:
                    break
                    
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    # Check for GFF3 version directive
                    if line.startswith('##gff-version 3'):
                        return Format.GFF3
                    continue
                
                # Check for GTF-specific attributes (gene_id and transcript_id)
                if re.search(r'gene_id\s+"[^"]+"\s*;\s*transcript_id\s+"[^"]+"', line):
                    return Format.GTF
                
                # Check for GFF3-specific attributes (ID=, Parent=)
                if re.search(r'ID=|Parent=', line):
                    return Format.GFF3
                
                # If we've reached a data line but couldn't determine format, assume GFF2
                return Format.GFF2
        
        # Default to GFF2 if we couldn't determine format
        return Format.GFF2