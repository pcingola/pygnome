"""
VCF Genotypes class for handling genotype data in VCF records.
"""
from dataclasses import dataclass
from typing import Any

from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType
from pygnome.parsers.vcf.vcf_format import VcfFormat
from pygnome.parsers.vcf.vcf_field_parser import decode_percent_encoded, encode_percent_encoded, parse_vcf_value


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


class VcfGenotypes:
    """
    Class for handling genotype data in VCF records.
    
    This class provides methods for accessing and modifying genotype data
    in VCF records. It delegates FORMAT field handling to VcfFormat.
    
    This class implements lazy parsing of genotype data to improve performance
    when only specific samples or fields are needed, especially for VCF files
    with millions of samples.
    """
    
    def __init__(self, format_and_genotypes_str: str, header: VcfHeader):
        """
        Initialize a VcfGenotypes object from a format and genotypes string.
        
        Args:
            format_and_genotypes_str: The tab-delimited string containing FORMAT and genotype fields
            header: The VCF header containing field definitions
        """
        self.header = header
        self._raw_str = format_and_genotypes_str
        self._fields: list[str] | None = None  # Will be populated on demand
        self._genotypes: list[Genotype] | None = None
        self._format: VcfFormat | None = None
        self._sample_data_cache: dict[int, dict[str, Any]] = {}  # Cache for parsed sample data
        self._modified_samples: dict[int, bool] = {}  # Track which samples have been modified
        
        # Check if we have genotype data (at least FORMAT field)
        self.has_genotypes = bool(format_and_genotypes_str.strip())
    
    def get_format(self) -> VcfFormat | None:
        """
        Get the FORMAT field object.
        
        Returns:
            The VcfFormat object, or None if there is no genotype data
        """
        if not self.has_genotypes:
            return None
        
        # Initialize format field if needed
        if self._format is None and self.has_genotypes:
            self._ensure_fields_parsed()
            format_str = self._fields[0] if self._fields else ""
            self._format = VcfFormat(format_str, self.header)
            
        return self._format

    def get_format_keys(self) -> list[str]:
        """
        Get the keys in the FORMAT field.
        
        Returns:
            A list of keys in the FORMAT field
        """
        if not self.has_genotypes:
            return []
            
        # Ensure format is initialized
        format_obj = self.get_format()
        if format_obj is None:
            return []
            
        return format_obj.get_keys()
    
    def has_format_key(self, key: str) -> bool:
        """
        Check if a key is present in the FORMAT field.
        
        Args:
            key: The key to check
            
        Returns:
            True if the key is present, False otherwise
        """
        if not self.has_genotypes:
            return False
            
        # Ensure format is initialized
        format_obj = self.get_format()
        if format_obj is None:
            return False
            
        return format_obj.has_key(key)
    
    def _ensure_fields_parsed(self) -> None:
        """
        Ensure the format and genotypes string has been parsed into fields.
        This implements lazy parsing - only parse when needed.
        """
        if self._fields is not None:
            return
            
        if not self._raw_str:
            self._fields = []
            return
            
        # Split the raw string into fields
        self._fields = self._raw_str.split("\t")
    
    def _parse_sample_data(self, sample_idx: int) -> dict[str, Any]:
        """
        Parse the data for a specific sample.
        
        Args:
            sample_idx: The index of the sample (0-based)
            
        Returns:
            A dictionary mapping FORMAT keys to values
        """
        # Check if we've already parsed this sample
        if sample_idx in self._sample_data_cache:
            return self._sample_data_cache[sample_idx]
        
        # Get the format keys
        format_keys = self.get_format_keys()
        if not format_keys:
            return {}
        
        # Ensure fields are parsed
        self._ensure_fields_parsed()
        
        # Check if we have enough fields for this sample
        if self._fields is None or sample_idx + 1 >= len(self._fields):
            return {}
        
        # Get the sample field (index 0 is FORMAT, so sample starts at index 1)
        sample_field = self._fields[sample_idx + 1]
        
        # Split the sample field into values
        sample_values = sample_field.split(":")
        
        # Create a dictionary mapping keys to values
        sample_data = {}
        for i, key in enumerate(format_keys):
            if i < len(sample_values):
                # Get the field definition
                field_def = self.header.get_format_field_definition(key)
                if field_def is None:
                    # Unknown field, store as string
                    sample_data[key] = sample_values[i]
                    continue
                
                # Parse the value based on the field type and number
                value = sample_values[i]
                
                if field_def is not None:
                    # Use the parse_vcf_value utility function
                    sample_data[key] = parse_vcf_value(value, field_def.type, field_def.number)
                else:
                    # Unknown field, store as string
                    sample_data[key] = value
            else:
                # Missing value
                sample_data[key] = None
        
        # Cache the parsed data
        self._sample_data_cache[sample_idx] = sample_data
        
        return sample_data
    
    def get_value(self, key: str, sample_idx: int) -> Any:
        """
        Get the value of a field for a specific sample.
        
        Args:
            key: The key of the field (must be in the FORMAT specification)
            sample_idx: The index of the sample (0-based)
            
        Returns:
            The value of the field, or None if not present
        """
        if not self.has_genotypes:
            return None
        
        # Check if the key is in the FORMAT field
        if not self.has_format_key(key):
            return None
        
        # Parse the sample data if needed
        sample_data = self._parse_sample_data(sample_idx)
        
        # Return the value
        return sample_data.get(key)
    
    def set_value(self, key: str, value: Any, sample_idx: int) -> None:
        """
        Set the value of a field for a specific sample.
        
        Args:
            key: The key of the field (must be in the FORMAT specification)
            value: The value to set
            sample_idx: The index of the sample (0-based)
        """
        if not self.has_genotypes:
            raise ValueError("Record has no genotype data")
        
        # Ensure fields are parsed
        self._ensure_fields_parsed()
        
        # Check if we have enough fields for the sample
        if self._fields is None or sample_idx + 1 >= len(self._fields):
            raise ValueError(f"Sample index out of range: {sample_idx}")
        
        # Check if the key is defined in the header
        field_def = self.header.get_format_field_definition(key)
        if field_def is None:
            raise ValueError(f"Unknown FORMAT field: {key}. Add it to the header first.")
        
        # Add the key to the FORMAT field if needed
        if not self.has_format_key(key) and self._format is not None:
            self._format.add_key(key)
        
        # Parse the sample data if needed
        if sample_idx not in self._sample_data_cache:
            self._parse_sample_data(sample_idx)
        
        # Update the cache
        if sample_idx in self._sample_data_cache:
            self._sample_data_cache[sample_idx][key] = value
        else:
            self._sample_data_cache[sample_idx] = {key: value}
        
        # Mark the sample as modified
        self._modified_samples[sample_idx] = True
        
        # Invalidate the raw string since we've made a modification
        self._raw_str = None
    
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
        
        # Check if GT is present in the FORMAT field
        if not self.has_format_key("GT"):
            # No genotype data
            self._genotypes = []
            return []
        
        # Parse the genotypes for all samples
        genotypes = []
        
        for i in range(len(self.header.samples)):
            # Get the GT value for this sample
            gt_value = self.get_value("GT", i)
            
            # If GT is missing or empty, add an empty genotype
            if gt_value is None or gt_value == "":
                genotypes.append(Genotype([], False))
                continue
            
            # Parse the genotype
            if "|" in gt_value:
                # Phased genotype
                allele_indices = []
                for allele in gt_value.split("|"):
                    if allele == ".":
                        allele_indices.append(None)
                    else:
                        try:
                            allele_indices.append(int(allele))
                        except ValueError:
                            allele_indices.append(None)
                genotypes.append(Genotype(allele_indices, True))
            else:
                # Unphased genotype
                allele_indices = []
                for allele in gt_value.split("/"):
                    if allele == ".":
                        allele_indices.append(None)
                    else:
                        try:
                            allele_indices.append(int(allele))
                        except ValueError:
                            allele_indices.append(None)
                genotypes.append(Genotype(allele_indices, False))
        
        # Cache the parsed genotypes
        self._genotypes = genotypes
        
        return genotypes
    
    def get_sample_names(self) -> list[str]:
        """
        Get the names of the samples in this record.
        """
        return self.header.samples
    
    def set_genotype(self, genotype: Genotype, sample_idx: int = 0) -> None:
        """
        Set the genotype for a specific sample.
        
        Args:
            genotype: The Genotype object
            sample_idx: The index of the sample (0-based)
        """
        if not self.has_genotypes:
            raise ValueError("Record has no genotype data")
        
        # Ensure fields are parsed
        self._ensure_fields_parsed()
        
        # Check if we have enough fields for the sample
        if self._fields is None or sample_idx + 1 >= len(self._fields):
            raise ValueError(f"Sample index out of range: {sample_idx}")
        
        # Update the genotype cache
        if self._genotypes is None:
            self._genotypes = self.get_genotypes()
        
        if self._genotypes and sample_idx < len(self._genotypes):
            self._genotypes[sample_idx] = genotype
        elif self._genotypes:
            # Extend the genotypes list if needed
            while len(self._genotypes) <= sample_idx:
                self._genotypes.append(Genotype([], False))
            self._genotypes[sample_idx] = genotype
        
        # Update the GT field in the sample data cache
        self.set_value("GT", str(genotype), sample_idx)
        
        # Note: set_value will invalidate the raw string
    
    def get_updated_format_and_genotypes(self) -> str:
        """
        Get the updated format and genotypes string.
        
        This method should be called before converting the record to a string.
        If no modifications have been made, returns the original raw string.
        
        Returns:
            The updated tab-delimited string containing FORMAT and genotype fields
        """
        if not self.has_genotypes:
            return ""
            
        # If no modifications have been made, return the original raw string
        if self._raw_str is not None:
            return self._raw_str
        
        # Ensure fields are parsed
        self._ensure_fields_parsed()
        
        if self._fields is None or not self._fields:
            return ""
        
        # Create a copy of the original fields
        updated_fields = self._fields.copy()
        
        # Update the FORMAT field
        if self._format is not None and self._format._modified:
            updated_fields[0] = self._format.to_string()
        
        # Update the sample fields
        for sample_idx, modified in self._modified_samples.items():
            if not modified:
                continue
            
            # Get the sample data
            sample_data = self._sample_data_cache.get(sample_idx, {})
            if not sample_data:
                continue
            
            # Get the format keys
            format_keys = self.get_format_keys()
            if not format_keys:
                continue
            
            # Build the sample field
            sample_values = []
            for key in format_keys:
                value = sample_data.get(key)
                if value is None:
                    sample_values.append(".")
                elif isinstance(value, str):
                    # Encode special characters
                    sample_values.append(encode_percent_encoded(value))
                elif isinstance(value, list):
                    # Format lists as comma-separated values
                    formatted_items = []
                    for item in value:
                        if item is None:
                            formatted_items.append(".")
                        else:
                            formatted_items.append(str(item))
                    sample_values.append(",".join(formatted_items))
                else:
                    sample_values.append(str(value))
            
            # Update the field
            field_idx = sample_idx + 1  # Index 0 is FORMAT, so samples start at index 1
            if field_idx < len(updated_fields):
                updated_fields[field_idx] = ":".join(sample_values)
            else:
                # Extend the fields array if needed
                while len(updated_fields) <= field_idx:
                    updated_fields.append(".")
                updated_fields[field_idx] = ":".join(sample_values)
        
        # Join the fields with tabs
        return "\t".join(updated_fields)
    
    def __str__(self) -> str:
        return self.get_updated_format_and_genotypes()
    
    def __repr__(self) -> str:
        return self.get_updated_format_and_genotypes()
