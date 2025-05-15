"""
MSI site record class for representing individual MSI site records.
"""
from dataclasses import dataclass


@dataclass
class MsiSiteRecord:
    """
    Represents a single record (line) in an MSI-sites file.
    
    Each record contains information about a microsatellite site in the genome,
    including its location, repeat unit, and flanking sequences.
    """
    
    chromosome: str
    location: int
    repeat_unit_length: int
    repeat_unit_binary: int
    repeat_times: int
    left_flank_binary: int
    right_flank_binary: int
    repeat_unit_bases: str
    left_flank_bases: str
    right_flank_bases: str
    
    def __str__(self) -> str:
        """Return a string representation of the MSI site record."""
        return (f"MsiSiteRecord(chromosome={self.chromosome}, location={self.location}, "
                f"repeat_unit={self.repeat_unit_bases}, repeat_times={self.repeat_times})")
    
    @property
    def end_location(self) -> int:
        """
        Calculate the end location of the microsatellite.
        
        Returns:
            The end position of the microsatellite (inclusive)
        """
        # The end location is the start location plus the length of the repeat
        # The length is the repeat unit length times the number of repeats
        return self.location + (self.repeat_unit_length * self.repeat_times) - 1
    
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