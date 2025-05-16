"""
VCF Record class for representing and parsing individual VCF records.
"""
import re
from dataclasses import dataclass
from typing import Any, Iterator
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand
from pygnome.genomics.variant import Variant
from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType, FieldNumber
from pygnome.parsers.vcf.variant_factory import VariantFactory


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


class VcfRecord(GenomicFeature):
    """
    Represents a single record (line) in a VCF file.
    
    This class implements lazy parsing of INFO and FORMAT fields to improve performance
    when only specific fields are needed.
    
    Inherits from GenomicFeature to provide standard genomic coordinate functionality.
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
        self._info_cache: dict[str, Any] = {}       # Cache for parsed INFO fields
        self._info_raw_cache: dict[str, str | None] | None = None  # Cache for raw string values
        self._info_modified: set[str] = set()  # Set of modified INFO fields
        self._format_cache: dict[str, list[Any]] = {}
        self._genotypes: list[Genotype] | None = None
        
        # Validate the number of fields
        if len(self._fields) < 8:
            raise ValueError(f"VCF record must have at least 8 fields: {line}")
        
        # Check if we have genotype data
        self.has_genotypes = len(self._fields) > 9
        
        # Initialize GenomicFeature fields
        chrom = self.get_chrom()
        start = self.get_pos()  # Already 0-based in get_pos()
        end = start + len(self.get_ref())  # End is exclusive
        id_value = self.get_id() if self.get_id() != "." else f"variant_{chrom}_{start}"
        
        # Call the parent class's __init__
        GenomicFeature.__init__(
            self,
            id=id_value,
            chrom=chrom,
            start=start,
            end=end,
            strand=Strand.UNSTRANDED
        )
    
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
    
    def _parse_info_raw(self, field_id: str) -> str | None:
        """Get the raw value of an INFO field."""
        if self._info_raw_cache is None:
            self._info_raw_cache = {}
            info_str = self._fields[7]
            if info_str != ".":
                # Split the INFO field into key-value pairs
                for field in info_str.split(";"):
                    if "=" in field:
                        key, value = field.split("=", 1)
                        self._info_raw_cache[key] = value
                    else:
                        # Flag field (presence means it's true)
                        self._info_raw_cache[field] = "True"
        # Return the raw value or None if not present
        return self._info_raw_cache.get(field_id)


    def get_info(self, field_id: str) -> Any:
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
        
        value_raw = self._parse_info_raw(field_id)
        if value_raw is None:
            # Field not present
            self._info_cache[field_id] = None
            return None
                
        # Get the field definition
        field_def = self.header.get_info_field_definition(field_id)
        if field_def is None:
            # Return None for unknown fields
            self._info_cache[field_id] = None
            return None
        
        # Handle flag fields
        if field_def.type == FieldType.FLAG:
            value_parsed = bool(value_raw)
        else:
            # Parse the value based on the field type and number
            value_parsed = self._parse_field_value(value_raw, field_def.type, field_def.number)
            
            # If we have a list with a single value, return the value directly for certain Number types
            if isinstance(value_parsed, list) and len(value_parsed) == 1:
                if field_def.number == 1 or field_def.number == FieldNumber.A and len(self.get_alt()) == 1:
                    value_parsed = value_parsed[0]
        
        # Cache the parsed value
        self._info_cache[field_id] = value_parsed
        
        return value_parsed
    
    def has_info(self, field_id: str) -> bool:
        """
        Check if an INFO field is present.
        """
        return self.get_info(field_id) is not None

    def get_format(self, field_id: str, sample_idx: int = 0) -> Any:
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
    
    def _parse_field_value(self, value: str, field_type: FieldType, number: FieldNumber | int | str) -> Any:
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
        # Removed debug print statement
        if number == 1:
            # Single value or unknown number
            return self._parse_single_value(value, field_type)
        elif number == 0:
            # Flag field (presence means it's true)
            return True
        else:
            # Multiple values - split by comma
            values = value.split(",")
            return [self._parse_single_value(v, field_type) for v in values]
    
    def _parse_single_value(self, value: str, field_type: FieldType) -> Any:
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
        
        if field_type == FieldType.INTEGER:
            return int(value)
        elif field_type == FieldType.FLOAT:
            return float(value)
        elif field_type == FieldType.CHARACTER:
            return value[0] if value else ""
        elif field_type == FieldType.FLAG:
            return True
        else:
            # String
            return decode_percent_encoded(value)
    
    def set_info(self, field_id: str, value: Any) -> None:
        """
        Set the value of an INFO field. If the field already exists, it will be updated.
        If it doesn't exist, it will be added.
        
        Args:
            field_id: The ID of the INFO field
            value: The value to set
        """
        # Mark the field as modified
        self._info_modified.add(field_id)
        # Get the field definition
        field_def = self.header.get_info_field_definition(field_id)
        if field_def is None:
            raise ValueError(f"Unknown INFO field: {field_id}. Add it to the header first.")
        # Split all raw info fields
        # Update the cache
        self._info_cache[field_id] = value
        # Invlidate the raw cache (first make sure the raw cache exists)
        self._parse_info_raw(field_id)
        self._info_raw_cache[field_id] = None  # type: ignore # Invalidate the raw cache
    
    def add_info(self, field_id: str, value: Any) -> None:
        """
        Add a new INFO field. If the field already exists, it will be updated.
        This is just an alias for set_info for clarity
        """
        self.set_info(field_id, value)
    
    def _info_str(self) -> str:
        """
        Create the string for an INFO field in the record.
        """
        # No modifications, return the original line
        if len(self._info_modified) == 0:
            return self._fields[7]

        # Re build the INFO field string
        assert self._info_raw_cache is not None
        infos = []
        for name, value in self._info_raw_cache.items():
            if name in self._info_modified:
                # Modified fields, use the new value, not the raw value
                value = self._info_cache[name]
                
                field_def = self.header.get_info_field_definition(name)
                if field_def is None:
                    raise ValueError(f"Unknown INFO field: {name}. Add it to the header first.")
                
                # Format the value based on its type
                formatted_value = self._format_info_value(value, field_def.type)
                infos.append(f"{name}={formatted_value}")
            else:
                # Unmodified fields, use the original 'raw' string'
                if value is not None:
                    infos.append(f"{name}={value}")
                else:
                    infos.append(f"{name}")

        return ";".join(infos)
    
    def _format_info_value(self, value: Any, field_type: FieldType) -> str:
        """
        Format a value for inclusion in an INFO field.
        
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
    
    def _update_raw_line(self) -> None:
        """Update the raw line to reflect changes to the fields."""
        self.raw_line = "\t".join(self._fields)
    
    # Methods for variant type detection
    def is_multi_allelic(self) -> bool:
        """
        Check if this variant has multiple alternate alleles.
        
        Returns:
            True if the variant has multiple alternate alleles, False otherwise
        """
        return len(self.get_alt()) > 1
    
    def get_end(self, alt_idx: int = 0) -> int:
        """
        Get the end position of the variant (0-based, exclusive).
        
        For SNPs, this is pos + 1. For indels, it's pos + len(ref).
        For structural variants, it's determined by the END or SVLEN INFO field.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            The end position (0-based, exclusive)
        """
        alt_alleles = self.get_alt()
        if alt_idx < len(alt_alleles) and VariantFactory.is_structural_variant(alt_alleles[alt_idx]):
            # Check for END info field
            if self.has_info("END"):
                end = self.get_info("END")
                if end is not None:
                    # END is 1-based inclusive in VCF, convert to 0-based exclusive
                    return end
            
            # Check for SVLEN info field
            if self.has_info("SVLEN"):
                svlen = self.get_info("SVLEN")
                if svlen is not None:
                    # SVLEN is the length of the variant
                    if isinstance(svlen, list):
                        # Use the SVLEN value for this alt allele if available
                        if alt_idx < len(svlen):
                            return self.get_pos() + abs(svlen[alt_idx])
                        # Otherwise use the first SVLEN value
                        return self.get_pos() + abs(svlen[0])
                    else:
                        return self.get_pos() + abs(svlen)
        
        # Default: pos + len(ref)
        return self.get_pos() + len(self.get_ref())
    
    def get_variant_str(self, alt_idx: int = 0) -> str:
        """
        Return a string representation of a specific variant.
        
        Args:
            alt_idx: Index of the alternate allele (default: 0)
            
        Returns:
            A string in the format "chrom:pos:ref>alt"
        """
        alt_alleles = self.get_alt()
        if alt_idx >= len(alt_alleles):
            raise IndexError(f"Alt index {alt_idx} out of range (0-{len(alt_alleles)-1})")
            
        return f"{self.get_chrom()}:{self.get_vcf_pos()}:{self.get_ref()}>{alt_alleles[alt_idx]}"
    
    def __iter__(self) -> Iterator[Variant]:
        """
        Iterate over the variants represented in this VCF record.
        
        For each alternate allele, yields a Variant object representing
        that specific variant.
        
        This allows iterating over all variants in a multi-allelic VCF record.
        
        Example:
            for variant in vcf_record:
                print(variant)
        """
        yield from VariantFactory.create_variants_from_record(self)

    def __str__(self) -> str:
        """Return the string representation of the VCF record."""
        if len(self._info_modified) == 0:
            # No modifications, return the original line
            return self.raw_line
        # If there are modifications, update the raw line

        
