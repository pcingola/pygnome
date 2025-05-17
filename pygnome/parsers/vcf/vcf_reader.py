"""
VCF Reader class for reading and parsing VCF files.
"""
import gzip
from pathlib import Path
from typing import Iterator, TextIO

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.vcf_field_parser import decode_percent_encoded


class VcfReader:
    """
    Reader for VCF (Variant Call Format) files.
    
    This class provides methods for reading and parsing VCF files, with support for
    both plain text and gzipped files. It implements lazy parsing of records to
    improve performance.
    
    The reader can be used as a context manager:
    
    ```python
    with VcfReader("example.vcf") as reader:
        for record in reader:
            # Process record
    ```
    
    Or manually:
    
    ```python
    reader = VcfReader("example.vcf")
    try:
        for record in reader:
            # Process record
    finally:
        reader.close()
    ```
    """
    
    def __init__(self, file_path: str | Path):
        """
        Initialize a VCF reader.
        
        Args:
            file_path: Path to the VCF file (can be plain text or gzipped)
        """
        self.file_path = Path(file_path)
        self.header = VcfHeader()
        self._file_handle: TextIO | None = None
        self._parse_header()
    
    def _open_file(self) -> TextIO:
        """
        Open the VCF file for reading.
        
        Returns:
            A file handle for the VCF file
        """
        if self.file_path.suffix == ".gz":
            return gzip.open(self.file_path, "rt")
        else:
            return open(self.file_path, "r")
    
    def _parse_header(self) -> None:
        """Parse the header section of the VCF file."""
        with self._open_file() as f:
            # Read lines until we find the header line
            for line in f:
                line = line.strip()
                
                if line.startswith("##"):
                    # Meta-information line
                    self.header.add_meta_line(line)
                elif line.startswith("#CHROM"):
                    # Header line
                    fields = line.split("\t")
                    
                    # Check if we have sample columns
                    if len(fields) > 8:
                        # Extract sample names (columns after FORMAT)
                        sample_names = fields[9:]
                        self.header.add_samples(sample_names)
                    
                    # We've reached the end of the header
                    break
    
    def __iter__(self) -> Iterator[VcfRecord]:
        """
        Iterate through the records in the VCF file.
        
        Yields:
            VcfRecord objects for each record in the file
        """
        # Close any existing file handle
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
        
        # Open the file
        self._file_handle = self._open_file()
        
        # Skip the header lines
        for line in self._file_handle:
            if not line.startswith("#"):
                # First data line
                yield VcfRecord(line, self.header)
                break
        
        # Read the rest of the file
        for line in self._file_handle:
            yield VcfRecord(line, self.header)
    
    def __enter__(self) -> "VcfReader":
        """
        Context manager entry.
        
        Returns:
            The VcfReader instance
        """
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Context manager exit.
        
        Ensures the file handle is closed when exiting the context.
        """
        self.close()
    
    def close(self) -> None:
        """Close the file handle."""
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
    
    def fetch(self, chrom: str, start: int, end: int) -> Iterator[VcfRecord]:
        """
        Fetch records in a specific genomic region.
        
        This is a simple implementation that iterates through all records and filters
        by position. For large files, an index should be used for efficient random access.
        
        Args:
            chrom: Chromosome name
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)
            
        Yields:
            VcfRecord objects for records in the specified region
        """
        for record in self:
            if record._parse_chrom() == chrom and record._parse_start() >= start and record._parse_start() < end:
                yield record
    
    def get_samples(self) -> list[str]:
        """
        Get the list of sample names.
        
        Returns:
            A list of sample names
        """
        return self.header.samples
    
    def get_contigs(self) -> list[str]:
        """
        Get the list of contigs.
        
        Returns:
            A list of contig names
        """
        return self.header.get_contigs()