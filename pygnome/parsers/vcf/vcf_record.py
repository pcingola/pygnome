"""
VCF Record class for representing and parsing individual VCF records.
"""
from dataclasses import dataclass
from pygnome.parsers.vcf.vcf_header import VcfHeader


@dataclass
class Genotype:
    """Represents a genotype from a VCF record."""
    allele_indices: list[int | None]
    phased: bool = False
    
    def __str__(self) -> str:
        """Return the string representation of the genotype."""
        separator = "|" if self.phased else "/"
        alleles = []
        for idx in self.allele_indices:
            if idx is None:
                alleles.append(".")
            else:
                alleles.append(str(idx))
        return separator.join(alleles)


class VcfRecord:
    """
    Represents a single record (line) in a VCF file.
    
    This class implements lazy parsing of INFO and FORMAT fields to improve performance
    when only specific fields are needed.
    """
    
    def __init__(self, line: str, header: VcfHeader):
        """
        Initialize a VCF record from a line in a VCF file.
        
        Args:
            line: A tab-delimited line from a VCF file
            header: The VCF header containing field definitions
        """
        self.raw_line = line
        self.header = header
        self._fields = line.strip().split("\t")
        
        # Cache for parsed fields
        self._info_cache: dict[str, any] = {}
        self._format_cache: dict[str, list[any]] = {}
        self._genotypes: list[Genotype] | None = None
        
        # Validate the number of fields
        if len(self._fields) < 8:
            raise ValueError(f"VCF record must have at least 8 fields: {line}")
        
        # Check if we have genotype data
        self.has_genotypes = len(self._fields) > 9
    
    def get_chrom(self) -> str:
        """Get the chromosome name."""
        return self._fields[0]
    
    def get_pos(self) -> int:
        """
        Get the position (0-based).
        
        VCF positions are 1-based, but we convert to 0-based internally.
        """
        return int(self._fields[1]) - 1
    
    def get_vcf_pos(self) -> int:
        """Get the original VCF position (1-based)."""
        return int(self._fields[1])
    
    def get_id(self) -> str:
        """Get the ID field."""
        return self._fields[2]
    
    def get_ref(self) -> str:
        """Get the reference allele."""
        return self._fields[3]
    
    def get_alt(self) -> list[str]:
        """Get the list of alternate alleles."""
        alt = self._fields[4]
        if alt == ".":
            return []
        return alt.split(",")
    
    def get_qual(self) -> float | None:
        """Get the quality score."""
        qual = self._fields[5]
        if qual == ".":
            return None
        return float(qual)
    
    def get_filter(self) -> list[str]:
        """Get the list of filters."""
        filter_field = self._fields[6]
        if filter_field == "." or filter_field == "PASS":
            return []
        return filter_field.split(";")
    
    def get_info(self, field_id: str) -> any:
        """
        Get the value of an INFO field.
        
        Args:
            field_id: The ID of the INFO field
            
        Returns:
            The parsed value of the INFO field, or None if not present
        """
        # Check if we've already parsed this field
        if field_id in self._info_cache:
            return self._info_cache[field_id]
        
        # Get the field definition
        field_def = self.header.get_info_field_definition(field_id)
        if field_def is None:
            # Return None for unknown fields
            self._info_cache[field_id] = None
            return None
        
        # Parse the INFO field
        info_str = self._fields[7]
        if info_str == ".":
            self._info_cache[field_id] = None
            return None
        
        # Split the INFO field into key-value pairs
        info_fields = {}
        for field in info_str.split(";"):
            if "=" in field:
                key, value = field.split("=", 1)
                info_fields[key] = value
            else:
                # Flag field (presence means it's true)
                info_fields[field] = True
        
        # Check if the field is present
        if field_id not in info_fields:
            self._info_cache[field_id] = None
            return None
        
        # Parse the field value based on its type
        value = info_fields[field_id]
        
        # Handle flag fields
        if field_def.type == "Flag":
            parsed_value = bool(value)
        else:
            # Parse the value based on the field type and number
            parsed_value = self._parse_field_value(value, field_def.type, field_def.number)
            
            # If we have a list with a single value, return the value directly for certain Number types
            if isinstance(parsed_value, list) and len(parsed_value) == 1:
                if field_def.number == "1" or field_def.number == "A" and len(self.get_alt()) == 1:
                    parsed_value = parsed_value[0]
        
        # Cache the parsed value
        self._info_cache[field_id] = parsed_value
        
        return parsed_value
    
    def has_info(self, field_id: str) -> bool:
        """
        Check if an INFO field is present.
        
        Args:
            field_id: The ID of the INFO field
            
        Returns:
            True if the field is present, False otherwise
        """
        info_str = self._fields[7]
        if info_str == ".":
            return False
        
        # Check if the field is present as a key=value pair or a flag
        for field in info_str.split(";"):
            if field == field_id or field.startswith(f"{field_id}="):
                return True
        
        return False
    
    def get_format(self, field_id: str, sample_idx: int = 0) -> any:
        """
        Get the value of a FORMAT field for a specific sample.
        
        Args:
            field_id: The ID of the FORMAT field
            sample_idx: The index of the sample (0-based)
            
        Returns:
            The parsed value of the FORMAT field, or None if not present
        """
        if not self.has_genotypes:
            return None
        
        # Check if we have enough fields for the sample
        if sample_idx >= len(self._fields) - 9:
            raise ValueError(f"Sample index out of range: {sample_idx}")
        
        # Check if we've already parsed this field for all samples
        if field_id in self._format_cache:
            return self._format_cache[field_id][sample_idx]
        
        # Get the field definition
        field_def = self.header.get_format_field_definition(field_id)
        if field_def is None:
            raise ValueError(f"Unknown FORMAT field: {field_id}")
        
        # Get the format keys and sample data
        format_keys = self._fields[8].split(":")
        if field_id not in format_keys:
            # Field not present in this record
            self._format_cache[field_id] = [None] * (len(self._fields) - 9)
            return None
        
        # Parse the field for all samples
        field_idx = format_keys.index(field_id)
        parsed_values = []
        
        for i in range(9, len(self._fields)):
            sample_data = self._fields[i].split(":")
            
            # Handle missing fields at the end
            if field_idx >= len(sample_data):
                parsed_values.append(None)
                continue
            
            # Get the raw value
            raw_value = sample_data[field_idx]
            if raw_value == ".":
                parsed_values.append(None)
                continue
            
            # Parse the value based on the field type and number
            parsed_value = self._parse_field_value(raw_value, field_def.type, field_def.number)
            parsed_values.append(parsed_value)
        
        # Cache the parsed values
        self._format_cache[field_id] = parsed_values
        
        return parsed_values[sample_idx]
    
    def has_format(self, field_id: str) -> bool:
        """
        Check if a FORMAT field is present.
        
        Args:
            field_id: The ID of the FORMAT field
            
        Returns:
            True if the field is present, False otherwise
        """
        if not self.has_genotypes:
            return False
        
        format_keys = self._fields[8].split(":")
        return field_id in format_keys
    
    def get_genotypes(self) -> list[Genotype]:
        """
        Get the genotypes for all samples.
        
        Returns:
            A list of Genotype objects, one for each sample
        """
        if not self.has_genotypes:
            return []
        
        # Check if we've already parsed the genotypes
        if self._genotypes is not None:
            return self._genotypes
        
        # Get the format keys and check if GT is present
        format_keys = self._fields[8].split(":")
        if "GT" not in format_keys:
            # No genotype data
            self._genotypes = []
            return []
        
        # Parse the genotypes for all samples
        gt_idx = format_keys.index("GT")
        genotypes = []
        
        for i in range(9, len(self._fields)):
            sample_data = self._fields[i].split(":")
            
            # Handle missing fields at the end
            if gt_idx >= len(sample_data):
                genotypes.append(Genotype([], False))
                continue
            
            # Get the raw genotype
            raw_gt = sample_data[gt_idx]
            if raw_gt == ".":
                genotypes.append(Genotype([], False))
                continue
            
            # Parse the genotype
            if "|" in raw_gt:
                # Phased genotype
                allele_indices = []
                for allele in raw_gt.split("|"):
                    if allele == ".":
                        allele_indices.append(None)
                    else:
                        allele_indices.append(int(allele))
                genotypes.append(Genotype(allele_indices, True))
            else:
                # Unphased genotype
                allele_indices = []
                for allele in raw_gt.split("/"):
                    if allele == ".":
                        allele_indices.append(None)
                    else:
                        allele_indices.append(int(allele))
                genotypes.append(Genotype(allele_indices, False))
        
        # Cache the parsed genotypes
        self._genotypes = genotypes
        
        return genotypes
    
    def get_sample_names(self) -> list[str]:
        """
        Get the names of the samples in this record.
        
        Returns:
            A list of sample names
        """
        return self.header.samples
    
    def get_alleles(self) -> list[str]:
        """
        Get all alleles (reference and alternates).
        
        Returns:
            A list of alleles, with the reference allele first
        """
        alleles = [self.get_ref()]
        alleles.extend(self.get_alt())
        return alleles
    
    def _parse_field_value(self, value: str, field_type: str, number: str) -> any:
        """
        Parse a field value based on its type and number.
        
        Args:
            value: The raw field value
            field_type: The type of the field (Integer, Float, String, Character, Flag)
            number: The number of values expected (1, A, R, G, ., etc.)
            
        Returns:
            The parsed field value
        """
        # Handle different number specifications
        if number == "1" or number == ".":
            # Single value or unknown number
            return self._parse_single_value(value, field_type)
        elif number == "0":
            # Flag field (presence means it's true)
            return True
        else:
            # Multiple values
            values = value.split(",")
            return [self._parse_single_value(v, field_type) for v in values]
    
    def _parse_single_value(self, value: str, field_type: str) -> any:
        """
        Parse a single field value based on its type.
        
        Args:
            value: The raw field value
            field_type: The type of the field (Integer, Float, String, Character, Flag)
            
        Returns:
            The parsed field value
        """
        if value == ".":
            return None
        
        if field_type == "Integer":
            return int(value)
        elif field_type == "Float":
            return float(value)
        elif field_type == "Character":
            return value[0] if value else ""
        else:
            # String or Flag
            return value
    
    def __str__(self) -> str:
        """Return the original VCF record line."""
        return self.raw_line