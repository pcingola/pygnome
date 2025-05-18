"""
ANN field parser for VCF files.

This module provides the AnnParser class for parsing the ANN field in VCF records,
which contains variant annotation information according to the VCF annotation format.
"""
from typing import Iterator, cast

from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.vcf_field_parser import decode_percent_encoded

from .variant_annotation import VariantAnnotation
from .enums import AnnotationImpact, FeatureType, ErrorWarningType, BiotypeCoding


class AnnParser:
    """
    Parser for the ANN field in VCF records.
    
    This class parses the ANN field in VCF records according to the VCF annotation
    format specification, yielding VariantAnnotation objects for each annotation.
    
    This class is both iterable and an iterator. To use it:
    
    1. Create an instance with a VCF record: parser = AnnParser(vcf_record)
    2. Iterate over the annotations: for annotation in parser: ...
    """
    
    # Field name in the VCF INFO field
    FIELD_NAME = "ANN"
    
    def __init__(self, record: VcfRecord):
        """
        Initialize the ANN parser with a VCF record.
        
        Args:
            record: A VcfRecord object containing the ANN field
        """
        self._record = record
        self._entries = []
        self._current_index = 0
        
    def parse(self):
        # Check if the record has the ANN field
        if self._entries:
            # Already parsed
            return
        if self._record and self._record.has_info(self.FIELD_NAME):
            # Get the ANN field value
            ann_value = self._record.get_info(self.FIELD_NAME)
            if ann_value:
                # The ANN field can contain multiple annotations
                if isinstance(ann_value, list):
                    # It's a list of entries
                    self._entries = ann_value
                else:
                    # It's a single entry
                    self._entries = [ann_value]
        else:
            self._entries = []
    
    def __iter__(self) -> "AnnParser":
        """
        Return self as an iterator.
        
        Returns:
            Self as an iterator
        """
        self._current_index = 0
        self.parse()
        return self
    
    def __next__(self) -> VariantAnnotation:
        """
        Get the next annotation.
        
        Returns:
            The next VariantAnnotation object
            
        Raises:
            StopIteration: When there are no more annotations
        """
        while self._current_index < len(self._entries):
            entry = self._entries[self._current_index]
            self._current_index += 1
            
            # Each annotation entry has fields separated by pipe characters
            fields = entry.split("|")
            
            # Skip if we don't have at least the required fields
            if len(fields) < 3:
                continue
            
            # Parse the fields
            annotation = self._parse_annotation_entry(fields)
            if annotation:
                return annotation
        
        raise StopIteration
    
    def _parse_annotation_entry(self, fields: list[str]) -> VariantAnnotation | None:
        """
        Parse a single annotation entry.
        
        Args:
            fields: A list of fields from a single annotation entry
            
        Returns:
            A VariantAnnotation object, or None if parsing failed
        """
        try:
            # Required fields
            allele = fields[0]
            annotation = fields[1]
            
            # Parse putative impact
            try:
                putative_impact = AnnotationImpact(fields[2])
            except ValueError:
                # Default to MODIFIER if the impact is invalid
                putative_impact = AnnotationImpact.MODIFIER
            
            # Create the annotation object with required fields
            result = VariantAnnotation(
                allele=allele,
                annotation=annotation,
                putative_impact=putative_impact
            )
            
            # Parse optional fields if available
            field_idx = 3
            
            # Gene Name
            if field_idx < len(fields) and fields[field_idx]:
                result.gene_name = decode_percent_encoded(fields[field_idx])
            field_idx += 1
            
            # Gene ID
            if field_idx < len(fields) and fields[field_idx]:
                result.gene_id = decode_percent_encoded(fields[field_idx])
            field_idx += 1
            
            # Feature Type
            if field_idx < len(fields) and fields[field_idx]:
                feature_type_str = decode_percent_encoded(fields[field_idx])
                try:
                    result.feature_type = FeatureType(feature_type_str)
                except ValueError:
                    # If not a standard feature type, use CUSTOM
                    result.feature_type = FeatureType.CUSTOM
            field_idx += 1
            
            # Feature ID
            if field_idx < len(fields) and fields[field_idx]:
                result.feature_id = decode_percent_encoded(fields[field_idx])
            field_idx += 1
            
            # Transcript Biotype
            if field_idx < len(fields) and fields[field_idx]:
                biotype_str = decode_percent_encoded(fields[field_idx])
                if biotype_str.lower() == "coding":
                    result.transcript_biotype = BiotypeCoding.CODING
                elif biotype_str.lower() == "noncoding":
                    result.transcript_biotype = BiotypeCoding.NONCODING
            field_idx += 1
            
            # Rank / Total
            if field_idx < len(fields) and fields[field_idx]:
                rank_total = fields[field_idx].split("/")
                if len(rank_total) == 2:
                    try:
                        result.rank = int(rank_total[0])
                        result.total = int(rank_total[1])
                    except ValueError:
                        pass
            field_idx += 1
            
            # HGVS.c
            if field_idx < len(fields) and fields[field_idx]:
                result.hgvs_c = decode_percent_encoded(fields[field_idx])
            field_idx += 1
            
            # HGVS.p
            if field_idx < len(fields) and fields[field_idx]:
                result.hgvs_p = decode_percent_encoded(fields[field_idx])
            field_idx += 1
            
            # cDNA_position / cDNA_len
            if field_idx < len(fields) and fields[field_idx]:
                cdna_pos_len = fields[field_idx].split("/")
                if len(cdna_pos_len) >= 1:
                    try:
                        result.cdna_pos = int(cdna_pos_len[0])
                        if len(cdna_pos_len) == 2:
                            result.cdna_length = int(cdna_pos_len[1])
                    except ValueError:
                        pass
            field_idx += 1
            
            # CDS_position / CDS_len
            if field_idx < len(fields) and fields[field_idx]:
                cds_pos_len = fields[field_idx].split("/")
                if len(cds_pos_len) >= 1:
                    try:
                        result.cds_pos = int(cds_pos_len[0])
                        if len(cds_pos_len) == 2:
                            result.cds_length = int(cds_pos_len[1])
                    except ValueError:
                        pass
            field_idx += 1
            
            # Protein_position / Protein_len
            if field_idx < len(fields) and fields[field_idx]:
                protein_pos_len = fields[field_idx].split("/")
                if len(protein_pos_len) >= 1:
                    try:
                        result.protein_pos = int(protein_pos_len[0])
                        if len(protein_pos_len) == 2:
                            result.protein_length = int(protein_pos_len[1])
                    except ValueError:
                        pass
            field_idx += 1
            
            # Distance to feature
            if field_idx < len(fields) and fields[field_idx]:
                try:
                    result.distance = int(fields[field_idx])
                except ValueError:
                    pass
            field_idx += 1
            
            # Errors, Warnings, or Information messages
            if field_idx < len(fields) and fields[field_idx]:
                message_codes = fields[field_idx].split("&")
                for code in message_codes:
                    message_type = ErrorWarningType.from_code(code.strip())
                    if message_type is not None:
                        result.add_message(message_type)
            
            return result
            
        except Exception as e:
            # If parsing fails, return None
            return None