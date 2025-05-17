"""
VCF Info class for handling INFO fields in VCF records.
"""
from typing import Any

from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType
from pygnome.parsers.vcf.vcf_field_parser import VcfFieldParser, encode_percent_encoded, decode_percent_encoded


class VcfInfo(VcfFieldParser):
    """
    Class for handling INFO fields in VCF records.
    
    This class provides methods for parsing, accessing, and modifying INFO fields
    in VCF records. It handles type conversion, percent encoding/decoding, and
    formatting of INFO fields according to the VCF specification.
    """
    
    def __init__(self, info_str: str, header: VcfHeader):
        """
        Initialize a VcfInfo object from an INFO field string.
        
        Args:
            info_str: The INFO field string from a VCF record
            header: The VCF header containing field definitions
        """
        super().__init__(info_str, header)
    
    def _ensure_parsed(self) -> None:
        """
        Ensure the INFO field string has been parsed.
        This implements lazy parsing - only parse when needed.
        """
        # If already parsed, return
        if self._field_raw_cache is not None:
            return
        
        # Initialize caches
        self._field_raw_cache = {}
        
        # Handle None or empty string
        if self._raw_str is None or self._raw_str == ".":
            return
        
        # Split the INFO field into key-value pairs
        for field in self._raw_str.split(";"):
            if "=" in field:
                key, value = field.split("=", 1)
                self._field_raw_cache[key] = value
            else:
                # Flag field (presence means it's true)
                self._field_raw_cache[field] = "True"
    
    def _get_field_definition(self, field_id: str) -> Any:
        """
        Get the definition for an INFO field.
        """
        return self.header.get_info_field_definition(field_id)
    
    def to_string(self) -> str:
        """
        Convert the INFO fields to a string representation.
        
        Returns:
            The string representation of the INFO fields
        """
        if not self.was_modified:
            # If the raw string is already set (i.e. not modified), return it
            return self._raw_str if self._raw_str is not None else "."

        self._ensure_parsed()

        # If no fields, return ".".
        if (self._field_raw_cache is None or not self._field_raw_cache):
            return "."

        # Build the INFO field string
        infos = []
        for name, value in self._field_raw_cache.items():
            if self.is_removed(name):
                # Skip deleted fields
                continue
            elif self.is_modified(name):
                # Modified fields, use the new value, not the raw value
                if self._field_cache is not None and name in self._field_cache:
                    value = self._field_cache[name]
                else:
                    # This should not happen, but just in case
                    raise ValueError(f"Field {name} not found in INFO cache.")
                
                field_def = self.header.get_info_field_definition(name)
                if field_def is None:
                    raise ValueError(f"Unknown INFO field: {name}. Add it to the header first.")
                
                # Skip flag fields with False value
                if field_def.type == FieldType.FLAG and not value:
                    continue
                
                # Format the value based on its type
                formatted_value = self._format_field_value(value, field_def.type)
                
                # Add the field to the list
                if field_def.type == FieldType.FLAG and value:
                    infos.append(name)
                else:
                    infos.append(f"{name}={formatted_value}")
            else:
                # Unmodified fields, use the original 'raw' string
                if value == "True":
                    infos.append(name)
                else:
                    infos.append(f"{name}={value}")
        
        # If no fields after processing, return "."
        if not infos:
            return "."
        
        return ";".join(infos)
    
    def __str__(self) -> str:
        """Return the string representation of the INFO fields."""
        return self.to_string()