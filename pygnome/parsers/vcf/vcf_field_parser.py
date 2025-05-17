"""
VCF Field Parser class for handling common parsing logic for VCF fields.
"""
import re
from typing import Any, Protocol

from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType, FieldNumber


def parse_vcf_value(value: str, field_type: FieldType, number: Any) -> Any:
    """
    Parse a VCF field value based on its type and number.
    
    Args:
        value: The raw field value
        field_type: The type of the field (Integer, Float, String, Character, Flag)
        number: The number of values expected (1, A, R, G, ., etc.)
        
    Returns:
        The parsed field value
    """
    # Handle missing value
    if value == ".":
        return None
    
    # Handle different number specifications
    if number == 1 or number == "1":
        # Single value
        return parse_single_value(value, field_type)
    elif number == 0 or number == "0":
        # Flag field (presence means it's true)
        return True
    elif "," in value:
        # Multiple values - split by comma
        values = value.split(",")
        return [parse_single_value(v, field_type) for v in values]
    else:
        # Single value but expected multiple - return as a list
        return [parse_single_value(value, field_type)]


def parse_single_value(value: str, field_type: FieldType) -> Any:
    """
    Parse a single VCF field value based on its type.
    
    Args:
        value: The raw field value
        field_type: The type of the field (Integer, Float, String, Character, Flag)
        
    Returns:
        The parsed field value
    """
    if value == ".":
        return None
    
    if field_type == FieldType.INTEGER:
        try:
            return int(value)
        except ValueError:
            return None
    elif field_type == FieldType.FLOAT:
        try:
            return float(value)
        except ValueError:
            return None
    elif field_type == FieldType.CHARACTER:
        return value[0] if value else ""
    elif field_type == FieldType.FLAG:
        return True
    else:
        # String
        return decode_percent_encoded(value)


def encode_percent_encoded(text: str) -> str:
    """
    Encode special characters in a string according to VCF specification.
    
    Args:
        text: The text to encode
        
    Returns:
        The encoded text with special characters percent-encoded
    """
    # We need to encode % first to avoid double-encoding
    text = text.replace('%', '%25')
    
    # Define a pattern to match special characters that need encoding
    special_chars = {
        ' ': '%20',  # Space
        ':': '%3A',
        ';': '%3B',
        '=': '%3D',
        ',': '%2C',
        '\r': '%0D',
        '\n': '%0A',
        '\t': '%09'
    }
    
    # Replace special characters with their percent-encoded equivalents
    result = text
    for char, encoded in special_chars.items():
        result = result.replace(char, encoded)
    
    return result


def decode_percent_encoded(text: str) -> str:
    """
    Decode percent-encoded characters in a string according to VCF specification.
    
    Args:
        text: The text to decode
        
    Returns:
        The decoded text
    """
    # Define a pattern to match percent-encoded characters
    pattern = re.compile(r'%([0-9A-Fa-f]{2})')
    
    # Define a replacement function
    def replace_match(match):
        hex_code = match.group(1)
        # Convert the hex code to an integer and then to a character
        return chr(int(hex_code, 16))
    
    # Replace all percent-encoded characters
    return pattern.sub(replace_match, text)


class FieldDefinitionProvider(Protocol):
    """Protocol for objects that can provide field definitions."""
    
    def get_info_field_definition(self, id: str) -> Any:
        """Get the definition for an INFO field."""
        ...
    
    def get_format_field_definition(self, id: str) -> Any:
        """Get the definition for a FORMAT field."""
        ...


class VcfFieldParser:
    """
    Base class for parsing VCF fields (INFO and FORMAT).
    
    This class provides common methods for parsing, accessing, and modifying
    fields in VCF records. It handles type conversion, percent encoding/decoding,
    and formatting of fields according to the VCF specification.
    
    This class implements lazy parsing - fields are only parsed when needed.
    """
    
    def __init__(self, raw_str: str, header: VcfHeader):
        """
        Initialize a VcfFieldParser object.
        
        Args:
            raw_str: The raw string containing the fields
            header: The VCF header containing field definitions
        """
        self.header = header
        self._raw_str = raw_str # Raw string representation of the fields (it is set to None if the fields are modified)
        self.was_modified: bool = False  # Flag to indicate if the field has been modified
        self._field_cache: dict[str, Any] | None = None  # Cache for parsed fields
        self._field_raw_cache: dict[str, str | None] | None = None  # Cache for raw string values
        self._field_modified: set[str] | None = None  # Set of modified fields
        self._field_removed: set[str] | None = None  # Set of deleted fields
    
    def _ensure_parsed(self) -> None:
        """
        Ensure the field string has been parsed.
        This implements lazy parsing - only parse when needed.
        
        This method should be implemented by derived classes to handle
        their specific field format.
        """
        raise NotImplementedError("Derived classes must implement _ensure_parsed")
    
    def _get_field_definition(self, field_id: str) -> Any:
        """
        Get the definition for a field.
        
        This method should be implemented by derived classes to return
        the appropriate field definition.
        
        Args:
            field_id: The ID of the field
            
        Returns:
            The field definition, or None if not found
        """
        raise NotImplementedError("Derived classes must implement _get_field_definition")
    
    
    def _format_field_value(self, value: Any, field_type: FieldType) -> str:
        """
        Format a value for inclusion in a field.
        
        Args:
            value: The value to format
            field_type: The type of the field
            
        Returns:
            The formatted value as a string
        """
        if value is None:
            return "."
        
        if isinstance(value, list):
            # Format each value in the list and join with commas
            formatted_values = []
            for v in value:
                if v is None:
                    formatted_values.append(".")
                elif field_type == FieldType.STRING:
                    formatted_values.append(encode_percent_encoded(str(v)))
                else:
                    formatted_values.append(str(v))
            return ",".join(formatted_values)
        else:
            # Format a single value
            if field_type == FieldType.STRING:
                return encode_percent_encoded(str(value))
            else:
                return str(value)
    
    def get(self, field_id: str) -> Any:
        """
        Get the value of a field.
        
        Args:
            field_id: The ID of the field
            
        Returns:
            The parsed value of the field, or None if not present
        """
        # Initialize field cache if needed
        if self._field_cache is None:
            self._field_cache = {}
        
        # Check if we've already parsed this field
        if field_id in self._field_cache:
            return self._field_cache[field_id]
        
        self._ensure_parsed()
        
        # If it's not parsed, check the raw cache and then parse it
        if self._field_raw_cache is None:
            return None
            
        value_raw = self._field_raw_cache.get(field_id)
        if value_raw is None:
            # Field not present
            self._field_cache[field_id] = None
            return None

        # Make sure the field has a definition
        field_def = self._get_field_definition(field_id)
        if field_def is None:
            # Return None for unknown fields
            self._field_cache[field_id] = None
            return None
        
        # Handle flag fields
        if field_def.type == FieldType.FLAG:
            value_parsed = bool(value_raw)
        else:
            # Parse the value based on the field type and number
            value_parsed = parse_vcf_value(value_raw, field_def.type, field_def.number)
            
            # If we have a list with a single value, return the value directly for certain Number types
            if isinstance(value_parsed, list) and len(value_parsed) == 1:
                if field_def.number == 1 or field_def.number == "1" or field_def.number == "A":
                    value_parsed = value_parsed[0]
        
        # Cache the parsed value
        self._field_cache[field_id] = value_parsed
        
        return value_parsed
    
    def has(self, field_id: str) -> bool:
        """
        Check if a field is present.
        """
        return self.get(field_id) is not None
    
    def remove(self, field_id: str) -> None:
        """
        Delete a field.
        """
        self._ensure_parsed()
        self._invalidate_cache(field_id, deleted=True)

    def _invalidate_cache(self, field_id: str, deleted: bool = False) -> None:
        """ Invalidate the raw string and mark as modified """
        self._raw_str = None    # Raw string representation of the fields (it is set to None if the fields are modified)
        self.was_modified = True  # Flag to indicate something was modified
        # Mark the field as deleted or modified
        if deleted:
            # Mark the field as deleted
            if self._field_removed is None:
                self._field_removed = set()
            self._field_removed.add(field_id)
        else:
            # Mark the field as modified
            if self._field_modified is None:
                self._field_modified = set()
            self._field_modified.add(field_id)
        # Invalidate the raw cache entry
        if self._field_raw_cache is not None:
            self._field_raw_cache[field_id] = None

    def is_removed(self, field_id: str) -> bool:
        """
        Check if a field has been removed.
        """
        return self._field_removed is not None and field_id in self._field_removed

    def is_modified(self, field_id: str) -> bool:
        """
        Check if a field has been modified.
        """
        return self._field_modified is not None and field_id in self._field_modified

    def set(self, field_id: str, value: Any) -> None:
        """
        Set the value of a field. If the field already exists, it will be updated.
        If it doesn't exist, it will be added.
        
        Args:
            field_id: The ID of the field
            value: The value to set
        """
        if value is None:
            # If the value is None, delete the field
            self.remove(field_id)
            return

        self._ensure_parsed()
        
        # Make sure the field has a definition
        field_def = self._get_field_definition(field_id)
        if field_def is None:
            raise ValueError(f"Unknown field: {field_id}. Add it to the header first.")
        
        # Update the field cache
        if self._field_cache is None:
            self._field_cache = {}    
        self._field_cache[field_id] = value
        
        # Invalidate caches
        self._invalidate_cache(field_id)
