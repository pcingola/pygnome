"""
VCF Record class for representing and parsing individual VCF records.
"""
from typing import Any, Iterator
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand
from pygnome.genomics.variant import Variant
from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType
from pygnome.parsers.vcf.variant_factory import VariantFactory
from pygnome.parsers.vcf.vcf_info import VcfInfo
from pygnome.parsers.vcf.vcf_format import VcfFormat
from pygnome.parsers.vcf.vcf_genotypes import VcfGenotypes, Genotype
from pygnome.parsers.vcf.vcf_field_parser import decode_percent_encoded


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
        
        # Split only the first 8 fields (fixed fields)
        parts = line.strip().split("\t", 8)  # Split into at most 9 parts
        
        # Validate the number of fields
        if len(parts) < 8:
            raise ValueError(f"VCF record must have at least 8 fields: {line}")
            
        # Store the first 8 fields
        self._fields = parts[:8]
        
        # Initialize GenomicFeature fields
        chrom = self._parse_chrom()
        start = self._parse_start()  # Already 0-based in get_pos()
        self.pos = start + 1  # Original VCF positions are 1-based
        end = start + len(self._parse_ref())  # End is exclusive
        id_value = self._parse_id()
        
        # Call the parent class's __init__
        GenomicFeature.__init__(self, id=id_value, chrom=chrom, start=start, end=end, strand=Strand.UNSTRANDED)

        self.ref = self._parse_ref()
        self.alts = self._parse_alt()
        self.qual = self._parse_qual()
        self.filter = self._parse_filter()

        # Initialize VcfInfo for handling INFO fields
        self.info = VcfInfo(self._fields[7], header)
        
        # Initialize VcfGenotypes for handling FORMAT and genotype fields
        if len(parts) > 8:
            # Pass the format and genotypes portion of the line (everything after INFO)
            format_and_genotypes_str = parts[8]
            self.genotypes = VcfGenotypes(format_and_genotypes_str, header)
        else:
            self.genotypes = VcfGenotypes("", header)
        

    def _parse_chrom(self) -> str:
        """Get the chromosome name."""
        return self._fields[0]
    
    def _parse_start(self) -> int:
        """
        Get the position (0-based).        
        VCF positions are 1-based, but we convert to 0-based internally.
        """
        return int(self._fields[1]) - 1
    
    def _parse_id(self) -> str:
        """Get the ID field."""
        _id = self._fields[2]
        if _id == ".":
            return ""
        return _id

    def _parse_ref(self) -> str:
        """Get the reference allele."""
        return self._fields[3]
    
    def _parse_alt(self) -> list[str]:
        """Get the list of alternate alleles."""
        alt = self._fields[4]
        if alt == ".":
            return []
        return alt.split(",")
    
    def _parse_qual(self) -> float | None:
        """Get the quality score."""
        qual = self._fields[5]
        if qual == ".":
            return None
        return float(qual)
    
    def _parse_filter(self) -> list[str]:
        """Get the list of filters."""
        filter_field = self._fields[6]
        if filter_field == "." or filter_field == "PASS":
            return []
        return filter_field.split(";")


    def get_info(self, field_id: str) -> Any:
        """
        Get the value of an INFO field.
        """
        return self.info.get(field_id)
    
    def has_info(self, field_id: str) -> bool:
        """
        Check if an INFO field is present.
        """
        return self.info.has(field_id)

    def get_format(self) -> VcfFormat | None:
        """
        Get the FORMAT field object.
        """
        return self.genotypes.get_format()

    def has_format(self, field_id: str) -> bool:
        """
        Check if a FORMAT field is present.
        """
        return self.genotypes.has_format_key(field_id)
    
    def get_genotypes(self) -> list[Genotype]:
        """
        Get the genotypes for all samples.
        """
        return self.genotypes.get_genotypes()
    
    def get_sample_names(self) -> list[str]:
        """
        Get the names of the samples in this record.
        """
        return self.genotypes.get_sample_names()
    
    def get_genotype_value(self, field_id: str, sample_idx: int = 0) -> Any:
        """
        Get the value of a genotype field for a specific sample.
        
        Args:
            field_id: The ID of the field (must be in the FORMAT specification)
            sample_idx: The index of the sample (0-based)
            
        Returns:
            The value of the field, or None if not present
        """
        return self.genotypes.get_value(field_id, sample_idx)
    
    def set_genotype_value(self, field_id: str, value: Any, sample_idx: int = 0) -> None:
        """
        Set the value of a genotype field for a specific sample.
        
        Args:
            field_id: The ID of the field (must be in the FORMAT specification)
            value: The value to set
            sample_idx: The index of the sample (0-based)
        """
        self.genotypes.set_value(field_id, value, sample_idx)
    
    def set_genotype(self, genotype: Genotype, sample_idx: int = 0) -> None:
        """
        Set the genotype for a specific sample.
        
        Args:
            genotype: The Genotype object
            sample_idx: The index of the sample (0-based)
        """
        self.genotypes.set_genotype(genotype, sample_idx)
    
    def get_alleles(self) -> list[str]:
        """
        Get all alleles (reference and alternates).
        
        Returns:
            A list of alleles, with the reference allele first
        """
        alleles = [self._parse_ref()]
        alleles.extend(self._parse_alt())
        return alleles
        
    def set_info(self, field_id: str, value: Any) -> None:
        """
        Set the value of an INFO field. If the field already exists, it will be updated.
        If it doesn't exist, it will be added.
        
        Args:
            field_id: The ID of the INFO field
            value: The value to set
        """
        self.info.set(field_id, value)
        
        # Invalidate the raw line
        self.raw_line = None
    
    def add_info(self, field_id: str, value: Any) -> None:
        """
        Add a new INFO field. If the field already exists, it will be updated.
        This is just an alias for set_info for clarity
        
        For flag fields, if value is False, the field is removed.
        """
        if value is False:
            field_def = self.header.get_info_field_definition(field_id)
            if field_def is not None and field_def.type == FieldType.FLAG:
                self.remove_info(field_id)
                return
        else:
            self.set_info(field_id, value)
    
    def remove_info(self, field_id: str) -> None:
        """
        Remove an INFO field.
        """
        self.info.remove(field_id)
    
    def _update_raw_line(self) -> None:
        """Update the raw line to reflect changes to the fields."""
        # Join the first 8 fields
        result = "\t".join(self._fields)
        
        # If we have genotype data, append the format and genotypes
        if self.genotypes.has_genotypes:
            # Get the updated format and genotypes string
            format_and_genotypes_str = self.genotypes.get_updated_format_and_genotypes()
            
            # If we have a format and genotypes string, append it to the result
            if format_and_genotypes_str:
                result += "\t" + format_and_genotypes_str
        
        self.raw_line = result
    
    # Methods for variant type detection
    def is_multi_allelic(self) -> bool:
        """
        Check if this variant has multiple alternate alleles.
        
        Returns:
            True if the variant has multiple alternate alleles, False otherwise
        """
        return len(self._parse_alt()) > 1
    
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
        alt_alleles = self._parse_alt()
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
                            return self._parse_start() + abs(svlen[alt_idx])
                        # Otherwise use the first SVLEN value
                        return self._parse_start() + abs(svlen[0])
                    else:
                        return self._parse_start() + abs(svlen)
        
        # Default: pos + len(ref)
        return self._parse_start() + len(self._parse_ref())
        
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
        # If no modifications have been made, return the original raw line
        if self.raw_line is not None and not self.info._field_modified and not self.genotypes._modified_samples:
            return self.raw_line
            
        # Update the INFO field
        self._fields[7] = self.info.to_string()
        
        # Join the first 8 fields
        result = "\t".join(self._fields)
        
        # If we have genotype data, append the format and genotypes
        if self.genotypes.has_genotypes:
            # Get the updated format and genotypes string
            format_and_genotypes_str = self.genotypes.get_updated_format_and_genotypes()
            
            # If we have a format and genotypes string, append it to the result
            if format_and_genotypes_str:
                result += "\t" + format_and_genotypes_str
        
        # Update the raw line
        self.raw_line = result
        
        return self.raw_line

        
