"""
Parser for FASTA format files.
"""

from dataclasses import dataclass
from pathlib import Path
import gzip
from typing import Iterator, Dict, Optional, TextIO

from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString


@dataclass
class FastaRecord:
    """
    Represents a single sequence record from a FASTA file.
    """
    identifier: str
    sequence: str
    description: str = ""
    
    def __str__(self) -> str:
        """Return a string representation of the record in FASTA format."""
        header = f">{self.identifier}"
        if self.description:
            header += f" {self.description}"
        
        # Format sequence with 80 characters per line
        formatted_seq = "\n".join(self.sequence[i:i+80] for i in range(0, len(self.sequence), 80))
        
        return f"{header}\n{formatted_seq}"


class FastaParser:
    """
    Parser for FASTA format files.
    """
    
    def __init__(self, file_path: Path):
        """Initialize the parser with a file path."""
        self.file_path = file_path
        self.file_handle = None
        self.current_header = None
        self.current_sequence = []
    
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
    
    def __next__(self) -> FastaRecord:
        """Get the next record from the file."""
        if not self.file_handle:
            raise RuntimeError("Parser not initialized. Use with statement.")
        
        while True:
            line = self.file_handle.readline()
            if not line:
                # End of file
                if self.current_header is not None:
                    # Return the last sequence
                    record = FastaRecord(
                        identifier=self.current_header.split()[0],
                        sequence=''.join(self.current_sequence),
                        description=' '.join(self.current_header.split()[1:]) if ' ' in self.current_header else ""
                    )
                    self.current_header = None
                    self.current_sequence = []
                    return record
                else:
                    # No more sequences
                    raise StopIteration
            
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            
            if line.startswith('>'):
                if self.current_header is not None:
                    # Return the previous sequence
                    record = FastaRecord(
                        identifier=self.current_header.split()[0],
                        sequence=''.join(self.current_sequence),
                        description=' '.join(self.current_header.split()[1:]) if ' ' in self.current_header else ""
                    )
                    self.current_header = line[1:]  # Remove the '>' prefix
                    self.current_sequence = []
                    return record
                else:
                    # Start a new sequence
                    self.current_header = line[1:]  # Remove the '>' prefix
                    self.current_sequence = []
            elif self.current_header is not None:
                # Add to the current sequence
                self.current_sequence.append(line)
    
    @staticmethod
    def parse(file_path: Path) -> Iterator[FastaRecord]:
        """
        Parse a FASTA file and yield FastaRecord objects.
        """
        with FastaParser(file_path) as parser:
            yield from parser
    
    @staticmethod
    def parse_as_dict(file_path: Path) -> Dict[str, str]:
        """Return a dictionary mapping identifiers to sequences."""
        return {record.identifier: record.sequence for record in FastaParser.parse(file_path)}
    
    @staticmethod
    def parse_first(file_path: Path) -> Optional[FastaRecord]:
        """Parse only the first sequence from a FASTA file."""
        try:
            return next(FastaParser.parse(file_path))
        except StopIteration:
            return None
    
    @staticmethod
    def parse_as_dna_strings(file_path: Path) -> Dict[str, DnaString]:
        """Return a dictionary mapping identifiers to DnaString objects."""
        return {record.identifier: DnaString(record.sequence) for record in FastaParser.parse(file_path)}
    
    @staticmethod
    def parse_as_rna_strings(file_path: Path) -> Dict[str, RnaString]:
        """Return a dictionary mapping identifiers to RnaString objects."""
        return {record.identifier: RnaString(record.sequence) for record in FastaParser.parse(file_path)}