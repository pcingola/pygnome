"""
Parser for FASTA format files.
"""

from dataclasses import dataclass
from pathlib import Path
import gzip
from typing import Iterator

from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString

DEFAULT_LENGTH_CONVERT_TO_DNA_STRING = 1024 * 1024  # 1 MB

@dataclass
class FastaRecord:
    """
    Represents a single sequence record from a FASTA file.
    """
    identifier: str
    sequence: str | DnaString
    description: str = ""
    
    def __str__(self) -> str:
        """Return a string representation of the record in FASTA format."""
        header = f">{self.identifier}"
        if self.description:
            header += f" {self.description}"
        
        # Format sequence with 80 characters per line
        formatted_seq = "\n".join(self.sequence[i:i+80] for i in range(0, len(self.sequence), 80))
        
        return f"{header}\n{formatted_seq}"
    
    @classmethod
    def create(cls, identifier: str, sequence: str, description: str = "", use_dna_string=False) -> 'FastaRecord':
        """
        Create a FastaRecord instance.
        """
        if use_dna_string:
            sequence = DnaString(sequence)
        return cls(identifier=identifier, sequence=sequence, description=description)


class FastaParser:
    """
    Parser for FASTA format files.
    """
    
    def __init__(self, file_path: Path, length_convert_to_dna_string: int = DEFAULT_LENGTH_CONVERT_TO_DNA_STRING):
        """Initialize the parser with a file path."""
        self.file_path = file_path
        self.file_handle = None
        self.current_header = None
        self.current_sequence = []        
        self.length_convert_to_dna_string = length_convert_to_dna_string

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
                    record = FastaRecord.create(
                        identifier=self.current_header.split()[0],
                        sequence=''.join(self.current_sequence),
                        description=' '.join(self.current_header.split()[1:]) if ' ' in self.current_header else "",
                        use_dna_string=len(self.current_sequence) > self.length_convert_to_dna_string
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
                    record = FastaRecord.create(
                        identifier=self.current_header.split()[0],
                        sequence=''.join(self.current_sequence),
                        description=' '.join(self.current_header.split()[1:]) if ' ' in self.current_header else "",
                        use_dna_string=len(self.current_sequence) > self.length_convert_to_dna_string
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
    
    def load(self) -> list[FastaRecord]:
        """
        Parse a FASTA file and yield FastaRecord objects.
        """
        records = []
        with self as parser:
            for record in parser:
                records.append(record)
        return records

    def load_as_dict(self) -> dict[str, FastaRecord]:
        """Return a dictionary mapping identifiers to FastaRecord objects."""
        result = {}
        with self as parser:
            for record in parser:
                result[record.identifier] = record
        return result
