"""
PyGnome parsers module for genomic annotation formats.

This module provides parsers for GFF2, GFF3, GTF, FASTA, and FASTQ formats with a unified
interface for working with genomic annotations and sequences.
"""

from .gff.record import GffRecord
from .gff.gff import Gff
from .gff.format import Format
from ..genomics.strand import Strand
from .gff.feature_hierarchy import FeatureHierarchy
from .fasta.fasta_parser import FastaParser, FastaRecord
from .fasta.fastq_parser import FastqParser, FastqRecord

__all__ = [
    'Record', 'Gff', 'Format', 'Strand', 'FeatureHierarchy',
    'FastaParser', 'FastaRecord', 'FastqParser', 'FastqRecord'
]