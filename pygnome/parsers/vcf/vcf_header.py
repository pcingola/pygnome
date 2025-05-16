"""
VCF Header class for parsing and representing VCF file headers.
"""
from dataclasses import dataclass, field
from enum import Enum, auto


class FieldType(str, Enum):
    """Enumeration of valid VCF field types."""
    INTEGER = "Integer"
    FLOAT = "Float"
    FLAG = "Flag"
    CHARACTER = "Character"
    STRING = "String"


class FieldNumber(str, Enum):
    """
    Enumeration of special VCF field number values.
    
    Values:
        A: One value per alternate allele
        R: One value per allele (including reference)
        G: One value per genotype
        DOT: Unknown or unbounded number of values
        LA: Local allele subset of A
        LR: Local allele subset of R
        LG: Local allele subset of G
        P: One value per allele in genotype
        M: One value per possible base modification
    """
    A = "A"  # One value per alternate allele
    R = "R"  # One value per allele (including reference)
    G = "G"  # One value per genotype
    DOT = "."  # Unknown or unbounded number of values
    LA = "LA"  # Local allele subset of A
    LR = "LR"  # Local allele subset of R
    LG = "LG"  # Local allele subset of G
    P = "P"  # One value per allele in genotype
    M = "M"  # One value per possible base modification
    
    @classmethod
    def from_str(cls, value: str) -> "FieldNumber | int | str":
        """
        Convert a string to a FieldNumber enum value or an integer.
        
        For numeric values, the integer value is returned.
        For special codes, the corresponding enum value is returned.
        
        Args:
            value: The string value to convert
            
        Returns:
            The corresponding FieldNumber enum value or an integer
        """
        try:
            # If it's a numeric value, return it as an integer
            return int(value)
        except ValueError:
            # If it's a special code, return the enum value
            try:
                return cls(value)
            except ValueError:
                # If it's not a valid enum value, return it as is
                # This should rarely happen with valid VCF files
                return value



@dataclass
class FieldDefinition:
    """Represents a structured field definition in the VCF header."""
    id: str
    number: FieldNumber | int | str
    type: FieldType
    description: str
    source: str | None = None
    version: str | None = None
    
    @classmethod
    def from_dict(cls, data: dict) -> "FieldDefinition":
        """Create a FieldDefinition from a dictionary of attributes."""
        # Convert number string to FieldNumber enum or keep as string for numeric values
        number_str = data.get("Number", "")
        number = FieldNumber.from_str(number_str)
        
        # Convert type string to FieldType enum
        type_str = data.get("Type", "")
        try:
            field_type = FieldType(type_str)
        except ValueError:
            # Default to STRING if the type is invalid
            field_type = FieldType.STRING
        
        return cls(
            id=data.get("ID", ""),
            number=number,
            type=field_type,
            description=data.get("Description", ""),
            source=data.get("Source"),
            version=data.get("Version")
        )


@dataclass
class FilterDefinition:
    """Represents a FILTER field definition in the VCF header."""
    id: str
    description: str
    
    @classmethod
    def from_dict(cls, data: dict) -> "FilterDefinition":
        """Create a FilterDefinition from a dictionary of attributes."""
        return cls(
            id=data.get("ID", ""),
            description=data.get("Description", "")
        )


@dataclass
class ContigDefinition:
    """Represents a contig definition in the VCF header."""
    id: str
    length: int | None = None
    md5: str | None = None
    url: str | None = None
    
    @classmethod
    def from_dict(cls, data: dict) -> "ContigDefinition":
        """Create a ContigDefinition from a dictionary of attributes."""
        length_str = data.get("length")
        length = int(length_str) if length_str is not None else None
        
        return cls(
            id=data.get("ID", ""),
            length=length,
            md5=data.get("md5"),
            url=data.get("URL")
        )


class VcfHeader:
    """
    Represents the header section of a VCF file.
    
    The header contains meta-information lines that describe the content and format
    of the VCF file, including INFO, FORMAT, FILTER, and contig definitions.
    """
    
    def __init__(self):
        """Initialize an empty VCF header."""
        self.meta_lines: list[str] = []
        self.file_format: str | None = None
        self.info_fields: dict[str, FieldDefinition] = {}
        self.format_fields: dict[str, FieldDefinition] = {}
        self.filters: dict[str, FilterDefinition] = {}
        self.contigs: dict[str, ContigDefinition] = {}
        self.alt_alleles: dict[str, str] = {}
        self.samples: list[str] = []
        self.other_metadata: dict[str, str] = {}
    
    def add_meta_line(self, line: str) -> None:
        """
        Add a meta-information line to the header.
        
        Args:
            line: A meta-information line starting with ##
        """
        if not line.startswith("##"):
            raise ValueError(f"Meta-information line must start with ##: {line}")
        
        self.meta_lines.append(line)
        self._parse_meta_line(line)
    
    def add_samples(self, sample_names: list[str]) -> None:
        """
        Add sample names from the header line.
        
        Args:
            sample_names: List of sample names
        """
        self.samples = sample_names
    
    def _parse_meta_line(self, line: str) -> None:
        """
        Parse a meta-information line and update the appropriate header section.
        
        Args:
            line: A meta-information line starting with ##
        """
        # Remove the ## prefix
        content = line[2:]
        
        # Check if it's a structured line (contains <...>)
        if "=<" in content and content.endswith(">"):
            key, value = content.split("=<", 1)
            value = value[:-1]  # Remove the trailing >
            self._parse_structured_meta(key, value)
        else:
            # It's an unstructured line
            if "=" in content:
                key, value = content.split("=", 1)
                if key == "fileformat":
                    self.file_format = value
                else:
                    self.other_metadata[key] = value
    
    def _parse_structured_meta(self, key: str, value: str) -> None:
        """
        Parse a structured meta-information line.
        
        Args:
            key: The key (e.g., INFO, FORMAT, FILTER)
            value: The value inside the angle brackets
        """
        # Parse the key-value pairs inside the angle brackets
        attributes = {}
        
        # Handle quoted values with commas
        parts = []
        current_part = ""
        in_quotes = False
        
        for char in value:
            if char == '"' and (not in_quotes or current_part[-1:] != "\\"):
                in_quotes = not in_quotes
                current_part += char
            elif char == "," and not in_quotes:
                parts.append(current_part)
                current_part = ""
            else:
                current_part += char
        
        if current_part:
            parts.append(current_part)
        
        # Parse each key=value pair
        for part in parts:
            if "=" in part:
                attr_key, attr_value = part.split("=", 1)
                # Remove quotes from quoted values
                if attr_value.startswith('"') and attr_value.endswith('"'):
                    attr_value = attr_value[1:-1]
                attributes[attr_key] = attr_value
        
        # Update the appropriate section based on the key
        if key == "INFO":
            field_def = FieldDefinition.from_dict(attributes)
            self.info_fields[field_def.id] = field_def
        elif key == "FORMAT":
            field_def = FieldDefinition.from_dict(attributes)
            self.format_fields[field_def.id] = field_def
        elif key == "FILTER":
            filter_def = FilterDefinition.from_dict(attributes)
            self.filters[filter_def.id] = filter_def
        elif key == "contig":
            contig_def = ContigDefinition.from_dict(attributes)
            self.contigs[contig_def.id] = contig_def
        elif key == "ALT":
            if "ID" in attributes and "Description" in attributes:
                self.alt_alleles[attributes["ID"]] = attributes["Description"]
    
    def get_info_field_definition(self, id: str) -> FieldDefinition | None:
        """
        Get the definition for an INFO field.
        
        Args:
            id: The ID of the INFO field
            
        Returns:
            The FieldDefinition for the INFO field, or None if not found
        """
        return self.info_fields.get(id)
    
    def get_format_field_definition(self, id: str) -> FieldDefinition | None:
        """
        Get the definition for a FORMAT field.
        
        Args:
            id: The ID of the FORMAT field
            
        Returns:
            The FieldDefinition for the FORMAT field, or None if not found
        """
        return self.format_fields.get(id)
    
    def get_contigs(self) -> list[str]:
        """
        Get the list of contigs defined in the header.
        
        Returns:
            A list of contig IDs
        """
        return list(self.contigs.keys())
    
    def __str__(self) -> str:
        """Return a string representation of the header."""
        lines = []
        
        # Add meta-information lines
        lines.extend(self.meta_lines)
        
        # Add the header line with sample names
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        if self.samples:
            header_line += "\tFORMAT\t" + "\t".join(self.samples)
        
        lines.append(header_line)
        
        return "\n".join(lines)