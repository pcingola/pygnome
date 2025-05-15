"""
VCF Record class for representing and parsing individual VCF records.
"""
from dataclasses import dataclass
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand
from pygnome.genomics.variant import VariantType, Variant, SNP, Insertion, Deletion, ComplexVariant
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
        self._info_cache: dict[str, any] = {}
        self._format_cache: dict[str, list[any]] = {}
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
    
    # Methods for variant type detection
    def is_snp(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is a SNP (Single Nucleotide Polymorphism).
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is a SNP, False otherwise
        """
        ref = self.get_ref()
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # SNP: single reference base and alternate allele is a single base
        return (len(ref) == 1 and len(alt) == 1 and alt in "ACGTN")
    
    def is_indel(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is an indel (insertion or deletion).
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is an indel, False otherwise
        """
        # Structural variants are not considered indels
        if self.is_structural_variant(alt_idx) or self.is_breakend(alt_idx):
            return False
            
        ref = self.get_ref()
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # Indel: reference and alternate allele have different lengths
        return len(ref) != len(alt)
    
    def is_insertion(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is an insertion.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is an insertion, False otherwise
        """
        # Structural variants are not considered insertions
        if self.is_structural_variant(alt_idx) or self.is_breakend(alt_idx):
            return False
            
        ref = self.get_ref()
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # Insertion: alternate allele is longer than the reference
        return len(alt) > len(ref)
    
    def is_deletion(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is a deletion.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is a deletion, False otherwise
        """
        # Structural variants are not considered deletions
        if self.is_structural_variant(alt_idx) or self.is_breakend(alt_idx):
            return False
            
        ref = self.get_ref()
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # Deletion: alternate allele is shorter than the reference
        return len(alt) < len(ref)
    
    def is_structural_variant(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is a structural variant.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is a structural variant, False otherwise
        """
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # Structural variant: alternate allele is symbolic (<...>)
        return alt.startswith("<") and alt.endswith(">")
    
    def is_breakend(self, alt_idx: int = 0) -> bool:
        """
        Check if this variant is a breakend.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            True if the variant is a breakend, False otherwise
        """
        alt_alleles = self.get_alt()
        
        if alt_idx >= len(alt_alleles):
            return False
            
        alt = alt_alleles[alt_idx]
        
        # Breakend: alternate allele contains [ or ]
        return "[" in alt or "]" in alt
    
    def is_multi_allelic(self) -> bool:
        """
        Check if this variant has multiple alternate alleles.
        
        Returns:
            True if the variant has multiple alternate alleles, False otherwise
        """
        return len(self.get_alt()) > 1
    
    def get_variant_type(self, alt_idx: int = 0) -> VariantType:
        """
        Get the type of this variant.
        
        Args:
            alt_idx: Index of the alternate allele to check (default: 0)
            
        Returns:
            A VariantType enum value describing the variant type
        """
        if self.is_snp(alt_idx):
            return VariantType.SNP
        elif self.is_insertion(alt_idx):
            return VariantType.INS
        elif self.is_deletion(alt_idx):
            return VariantType.DEL
        elif self.is_structural_variant(alt_idx):
            # Get the specific type from the symbolic allele
            alt_alleles = self.get_alt()
            if alt_idx < len(alt_alleles):
                alt = alt_alleles[alt_idx]
                if alt.startswith("<") and alt.endswith(">"):
                    # Try to map the symbolic allele to a variant type
                    sv_type = alt[1:-1]  # Remove the < and >
                    return sv_type
            return VariantType.SV
        elif self.is_breakend(alt_idx):
            return VariantType.BND
        else:
            return VariantType.OTHER
    
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
        if self.is_structural_variant(alt_idx):
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
    
    def __iter__(self):
        """
        Iterate over the variants represented in this VCF record.
        
        For each alternate allele, yields a Variant object representing
        that specific variant.
        
        This allows iterating over all variants in a multi-allelic VCF record.
        
        Example:
            for variant in vcf_record:
                print(variant)
        """
        for i in range(len(self.get_alt())):
            variant_id = self.get_id() if self.get_id() != "." else f"variant_{self.get_chrom()}_{self.get_pos()}"
            info = {k: self.get_info(k) for k in self.header.info_fields if self.has_info(k)}
            
            if self.is_snp(i):
                yield SNP(id=variant_id, chrom=self.get_chrom(), start=self.get_pos(), end=self.get_end(i), strand=Strand.UNSTRANDED, ref=self.get_ref(), alt=self.get_alt()[i], quality=self.get_qual(), filters=self.get_filter(), info=info)
            elif self.is_insertion(i):
                yield Insertion(id=variant_id, chrom=self.get_chrom(), start=self.get_pos(), end=self.get_end(i), strand=Strand.UNSTRANDED, ref=self.get_ref(), alt=self.get_alt()[i], quality=self.get_qual(), filters=self.get_filter(), info=info)
            elif self.is_deletion(i):
                yield Deletion(id=variant_id, chrom=self.get_chrom(), start=self.get_pos(), end=self.get_end(i), strand=Strand.UNSTRANDED, ref=self.get_ref(), alt=self.get_alt()[i], quality=self.get_qual(), filters=self.get_filter(), info=info)
            else:
                yield ComplexVariant(id=variant_id, chrom=self.get_chrom(), start=self.get_pos(), end=self.get_end(i), strand=Strand.UNSTRANDED, ref=self.get_ref(), alt=self.get_alt()[i], quality=self.get_qual(), filters=self.get_filter(), info=info, description=f"{self.get_variant_type(i)} variant")
