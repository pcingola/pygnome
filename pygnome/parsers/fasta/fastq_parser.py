"""
Parser for FASTQ format files.
"""

from dataclasses import dataclass
from pathlib import Path
import gzip
from typing import Iterator, Dict, Optional, List, Tuple, TextIO

from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString


@dataclass
class FastqRecord:
    """
    Represents a single sequence record from a FASTQ file.
    """
    identifier: str
    sequence: str
    quality: str
    description: str = ""
    
    def __post_init__(self):
        """Validate that sequence and quality have the same length."""
        if len(self.sequence) != len(self.quality):
            raise ValueError(f"Sequence and quality must have the same length: {len(self.sequence)} != {len(self.quality)}")
    
    def __str__(self) -> str:
        """Return a string representation of the record in FASTQ format."""
        header = f"@{self.identifier}"
        if self.description:
            header += f" {self.description}"
        
        return f"{header}\n{self.sequence}\n+\n{self.quality}"
    
    def get_quality_scores(self) -> List[int]:
        """Convert ASCII quality characters to Phred quality scores."""
        return [ord(c) - 33 for c in self.quality]  # Assuming Sanger format (Phred+33)


class FastqParser:
    """
    Parser for FASTQ format files.
    """
    
    def __init__(self, file_path: Path):
        """Initialize the parser with a file path."""
        self.file_path = file_path
        self.file_handle = None
    
    def __enter__(self):
        """Context manager entry point."""
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {self.file_path}")
        
        # Open file (handle gzipped files)
        open_func = gzip.open if str(self.file_path).endswith(('.gz', '.gzip')) else open
        self.file_handle = open_func(self.file_path, 'rt')
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point."""
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None
    
    def __iter__(self):
        """Return self as iterator."""
        return self
    
    def __next__(self) -> FastqRecord:
        """Get the next record from the file."""
        if not self.file_handle:
            raise RuntimeError("Parser not initialized. Use with statement.")
        
        # Read the four lines that make up a FASTQ record
        header = self.file_handle.readline().strip()
        if not header:
            raise StopIteration  # End of file
        
        if not header.startswith('@'):
            raise ValueError(f"Invalid FASTQ format: header line does not start with '@': {header}")
        
        sequence = self.file_handle.readline().strip()
        if not sequence:
            raise ValueError("Invalid FASTQ format: missing sequence line")
        
        plus_line = self.file_handle.readline().strip()
        if not plus_line.startswith('+'):
            raise ValueError(f"Invalid FASTQ format: third line does not start with '+': {plus_line}")
        
        quality = self.file_handle.readline().strip()
        if not quality:
            raise ValueError("Invalid FASTQ format: missing quality line")
        
        # Parse the header
        header = header[1:]  # Remove the '@' prefix
        parts = header.split(maxsplit=1)
        identifier = parts[0]
        description = parts[1] if len(parts) > 1 else ""
        
        return FastqRecord(identifier, sequence, quality, description)
    
    @staticmethod
    def parse(file_path: Path) -> Iterator[FastqRecord]:
        """Parse a FASTQ file and yield FastqRecord objects."""
        with FastqParser(file_path) as parser:
            yield from parser
    
    @staticmethod
    def parse_as_dict(file_path: Path) -> Dict[str, Tuple[str, str]]:
        """Return a dictionary mapping identifiers to (sequence, quality) tuples."""
        return {record.identifier: (record.sequence, record.quality) 
                for record in FastqParser.parse(file_path)}
    
    @staticmethod
    def parse_first(file_path: Path) -> Optional[FastqRecord]:
        """Parse only the first sequence from a FASTQ file."""
        try:
            return next(FastqParser.parse(file_path))
        except StopIteration:
            return None
    
    @staticmethod
    def parse_as_dna_strings(file_path: Path) -> Dict[str, Tuple[DnaString, str]]:
        """Return a dictionary mapping identifiers to (DnaString, quality) tuples."""
        return {record.identifier: (DnaString(record.sequence), record.quality) 
                for record in FastqParser.parse(file_path)}
    
    @staticmethod
    def parse_as_rna_strings(file_path: Path) -> Dict[str, Tuple[RnaString, str]]:
        """Return a dictionary mapping identifiers to (RnaString, quality) tuples."""
        return {record.identifier: (RnaString(record.sequence), record.quality) 
                for record in FastqParser.parse(file_path)}