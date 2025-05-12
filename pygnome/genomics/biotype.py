"""Biotype enumeration for genomic features."""

from enum import Enum


class Biotype(str, Enum):
    """Common biotypes for genes and transcripts."""
    PROTEIN_CODING = "protein_coding"
    LNCRNA = "lncRNA"
    MIRNA = "miRNA"
    SNRNA = "snRNA"
    SNORNA = "snoRNA"
    RRNA = "rRNA"
    TRNA = "tRNA"
    PSEUDOGENE = "pseudogene"
    PROCESSED_PSEUDOGENE = "processed_pseudogene"
    UNPROCESSED_PSEUDOGENE = "unprocessed_pseudogene"
    NONSENSE_MEDIATED_DECAY = "nonsense_mediated_decay"
    RETAINED_INTRON = "retained_intron"
    PROCESSED_TRANSCRIPT = "processed_transcript"
    UNKNOWN = "unknown"