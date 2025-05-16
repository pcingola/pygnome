"""
MSI sites reader class for reading and parsing MSI sites files.
"""
import gzip
from pathlib import Path
from typing import Iterator, TextIO, Optional

from .msi_site_record import MsiSiteRecord


class MsiSitesReader:
    """
    Reader for MSI (Microsatellite Instability) sites files.
    
    This class provides methods for reading and parsing MSI sites files, with support for
    both plain text and gzipped files. It implements lazy parsing of records to
    improve performance when dealing with large files (e.g., 2GB files with 30+ million lines).
    
    The reader can be used as a context manager:
    
    ```python
    with MsiSitesReader("example.txt") as reader:
        for record in reader:
            # Process record
    ```
    
    Or manually:
    
    ```python
    reader = MsiSitesReader("example.txt")
    try:
        for record in reader:
            # Process record
    finally:
        reader.close()
    ```
    """
    
    def __init__(self, file_path: Path):
        """
        Initialize an MSI sites reader.
        
        Args:
            file_path: Path to the MSI sites file (can be plain text or gzipped)
        """
        self.file_path = file_path
        self._file_handle: TextIO | None = None
        self._header_processed = False
        self._column_indices = {}
        
        # Check if file exists
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
    
    def _open_file(self) -> TextIO:
        """
        Open the MSI sites file for reading.
        
        Returns:
            A file handle for the MSI sites file
        """
        if str(self.file_path).endswith(('.gz', '.gzip')):
            return gzip.open(self.file_path, "rt")
        else:
            return open(self.file_path, "r")
    
    def _process_header(self, header_line: str) -> None:
        """
        Process the header line to determine column indices.
        
        Args:
            header_line: The header line from the MSI sites file
        """
        columns = header_line.strip().split('\t')
        
        # Map column names to indices
        for i, column in enumerate(columns):
            self._column_indices[column] = i
        
        # Verify that all required columns are present
        required_columns = [
            "chromosome", "location", "repeat_unit_length", "repeat_unit_binary",
            "repeat_times", "left_flank_binary", "right_flank_binary",
            "repeat_unit_bases", "left_flank_bases", "right_flank_bases"
        ]
        
        missing_columns = [col for col in required_columns if col not in self._column_indices]
        if missing_columns:
            raise ValueError(f"Missing required columns in MSI sites file: {', '.join(missing_columns)}")
        
        self._header_processed = True
    
    def _parse_line(self, line: str) -> MsiSiteRecord:
        """
        Parse a single line from the MSI sites file.
        
        Args:
            line: A tab-delimited line from the MSI sites file
            
        Returns:
            An MsiSiteRecord object representing the parsed line
        """
        fields = line.strip().split('\t')
        
        # Ensure we have enough fields
        if len(fields) < len(self._column_indices):
            raise ValueError(f"Line has fewer fields than expected: {line}")
        
        # Extract fields using column indices
        return MsiSiteRecord(
            chrom=fields[self._column_indices["chromosome"]],  # Use chrom instead of chromosome
            location=int(fields[self._column_indices["location"]]),
            repeat_unit_length=int(fields[self._column_indices["repeat_unit_length"]]),
            repeat_unit_binary=int(fields[self._column_indices["repeat_unit_binary"]]),
            repeat_times=int(fields[self._column_indices["repeat_times"]]),
            left_flank_binary=int(fields[self._column_indices["left_flank_binary"]]),
            right_flank_binary=int(fields[self._column_indices["right_flank_binary"]]),
            repeat_unit_bases=fields[self._column_indices["repeat_unit_bases"]],
            left_flank_bases=fields[self._column_indices["left_flank_bases"]],
            right_flank_bases=fields[self._column_indices["right_flank_bases"]]
        )
    
    def __iter__(self) -> Iterator[MsiSiteRecord]:
        """
        Iterate through the records in the MSI sites file.
        
        Yields:
            MsiSiteRecord objects for each record in the file
        """
        # Close any existing file handle
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
        
        # Open the file
        self._file_handle = self._open_file()
        
        # Process the header if not already done
        if not self._header_processed:
            header_line = next(self._file_handle, None)
            if header_line is None:
                raise ValueError("Empty MSI sites file")
            self._process_header(header_line)
        else:
            # Skip the header line
            next(self._file_handle, None)
        
        # Read and parse each line
        for line in self._file_handle:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            
            try:
                yield self._parse_line(line)
            except Exception as e:
                # Log error but continue parsing
                print(f"Error parsing line: {line}")
                print(f"Error: {str(e)}")
    
    def __enter__(self) -> "MsiSitesReader":
        """
        Context manager entry.
        
        Returns:
            The MsiSitesReader instance
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
    
    def fetch(self, chromosome: str, start: int, end: int) -> Iterator[MsiSiteRecord]:
        """
        Fetch records in a specific genomic region.
        
        This is a simple implementation that iterates through all records and filters
        by position. For large files, an index should be used for efficient random access.
        
        Args:
            chromosome: Chromosome name
            start: Start position (0-based, inclusive)
            end: End position (0-based, inclusive)
            
        Yields:
            MsiSiteRecord objects for records in the specified region
        """
        for record in self:
            # Check if the record overlaps with the specified region
            if (record.chrom == chromosome and
                record.location <= end and
                record.end_location >= start):
                yield record