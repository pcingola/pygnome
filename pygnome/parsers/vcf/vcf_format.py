"""
VCF Format class for handling FORMAT fields in VCF records.
"""
from typing import List

from pygnome.parsers.vcf.vcf_header import VcfHeader


class VcfFormat:
    """
    Class for handling FORMAT fields in VCF records.
    
    This class provides methods for accessing and modifying the FORMAT field
    in VCF records. The FORMAT field specifies the format of the genotype fields
    for all samples but doesn't contain values itself.
    """
    
    def __init__(self, format_str: str, header: VcfHeader):
        """
        Initialize a VcfFormat object from a FORMAT field string.
        
        Args:
            format_str: The FORMAT field string from a VCF record (e.g., "GT:GQ:DP:HQ")
            header: The VCF header containing field definitions
        """
        self.header = header
        self._raw_str = format_str
        self._keys: List[str] | None = None
        self._modified = False
    
    def _ensure_parsed(self) -> None:
        """
        Ensure the FORMAT field string has been parsed.
        This implements lazy parsing - only parse when needed.
        """
        # If already parsed, return
        if self._keys is not None:
            return
        
        # Handle None or empty string
        if self._raw_str is None or self._raw_str == ".":
            self._keys = []
            return
        
        # Split the FORMAT field into keys
        self._keys = self._raw_str.split(":")
    
    def get_keys(self) -> List[str]:
        """
        Get the keys in the FORMAT field.
        
        Returns:
            A list of keys in the FORMAT field
        """
        self._ensure_parsed()
        return self._keys.copy() if self._keys else []
    
    def has_key(self, key: str) -> bool:
        """
        Check if a key is present in the FORMAT field.
        
        Args:
            key: The key to check
            
        Returns:
            True if the key is present, False otherwise
        """
        self._ensure_parsed()
        return key in self.get_keys()
    
    def add_key(self, key: str) -> None:
        """
        Add a key to the FORMAT field.
        
        Args:
            key: The key to add
        """
        self._ensure_parsed()
        
        # Check if the key already exists
        if key in self.get_keys():
            return
        
        # Check if the key is defined in the header
        field_def = self.header.get_format_field_definition(key)
        if field_def is None:
            raise ValueError(f"Unknown FORMAT field: {key}. Add it to the header first.")
        
        # Add the key
        if self._keys is None:
            self._keys = []
        self._keys.append(key)
        
        # Mark as modified
        self._modified = True
    
    def remove_key(self, key: str) -> None:
        """
        Remove a key from the FORMAT field.
        
        Args:
            key: The key to remove
        """
        self._ensure_parsed()
        
        # Check if the key exists
        if self._keys is None or key not in self._keys:
            return
        
        # Remove the key
        self._keys.remove(key)
        
        # Mark as modified
        self._modified = True
    
    def to_string(self) -> str:
        """
        Convert the FORMAT field to a string representation.
        
        Returns:
            The string representation of the FORMAT field
        """
        self._ensure_parsed()
        
        # If no keys or empty string, return "."
        if not self._keys or self._raw_str == "":
            return "."
        
        # Join the keys with colons
        return ":".join(self._keys)
    
    def __str__(self) -> str:
        """Return the string representation of the FORMAT field."""
        return self.to_string()