"""
MSI site record class for representing individual MSI site records.
"""
from dataclasses import dataclass, field

from ...genomics.genomic_feature import GenomicFeature
from ...genomics.strand import Strand


class MsiSiteRecord(GenomicFeature):
    """
    Represents a single record (line) in an MSI-sites file.
    
    Each record contains information about a microsatellite site in the genome,
    including its location, repeat unit, and flanking sequences.
    
    Inherits from GenomicFeature to provide standard genomic coordinate functionality.
    """
    
    def __init__(
        self,
        chrom: str,
        location: int,
        repeat_unit_length: int,
        repeat_times: int,
        repeat_unit_bases: str,
        id: str = "",
        left_flank_bases: str = '',
        right_flank_bases: str = '',
        left_flank_binary: int = 0,
        right_flank_binary: int = 0,
        repeat_unit_binary: int = 0,
        strand: Strand = Strand.UNSTRANDED
    ):
        """
        Initialize an MSI site record.
        
        Args:
            repeat_unit_length: Length of the repeat unit
            repeat_unit_binary: Binary representation of the repeat unit
            repeat_times: Number of times the repeat unit is repeated
            left_flank_binary: Binary representation of the left flanking sequence
            right_flank_binary: Binary representation of the right flanking sequence
            repeat_unit_bases: Base sequence of the repeat unit
            left_flank_bases: Base sequence of the left flanking region
            right_flank_bases: Base sequence of the right flanking region
            chromosome: Chromosome name
            location: Start location of the microsatellite
            id: Optional identifier for the feature
            strand: Strand of the feature
        """
        # Store MSI-specific fields
        self.repeat_unit_length = repeat_unit_length
        self.repeat_unit_binary = repeat_unit_binary
        self.repeat_times = repeat_times
        self.left_flank_binary = left_flank_binary
        self.right_flank_binary = right_flank_binary
        self.repeat_unit_bases = repeat_unit_bases
        self.left_flank_bases = left_flank_bases
        self.right_flank_bases = right_flank_bases
        self.location = location
        
        # Calculate GenomicFeature fields
        start = location
        end = location + (repeat_unit_length * repeat_times)
        
        # Set a default ID if not provided
        if not id:
            id = f"MSI_{chrom}_{start}"
        
        # Initialize the parent class
        super().__init__(id=id, chrom=chrom, start=start, end=end, strand=strand)
    
    # __post_init__ is not needed as we handle initialization in __init__
    
    @property
    def end_location(self) -> int:
        """
        Calculate the end location of the microsatellite.
        
        Returns:
            The end position of the microsatellite (inclusive)
        """
        # The end location is the start location plus the length of the repeat
        # The length is the repeat unit length times the number of repeats
        # Note: This returns the inclusive end position for backward compatibility
        return self.end - 1
    
    @property
    def sequence(self) -> str:
        """
        Get the complete sequence of the microsatellite.
        
        Returns:
            The microsatellite sequence (repeat unit repeated n times)
        """
        return self.repeat_unit_bases * self.repeat_times
    
    @property
    def full_sequence(self) -> str:
        """
        Get the full sequence including flanking regions.
        
        Returns:
            The full sequence: left flank + microsatellite + right flank
        """
        return self.left_flank_bases + self.sequence + self.right_flank_bases