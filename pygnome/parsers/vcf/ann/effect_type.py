"""
Effect Type enumeration for VCF annotations.

This module defines the EffectType enumeration used in parsing the ANN field in VCF files,
representing the specific effects a variant can have on genomic features.
"""
from enum import Enum
import warnings

from .enums import AnnotationImpact


class EffectType(str, Enum):
    """
    Enumeration of effect types for variant annotations.
    
    These types represent the specific effects a variant can have on genomic features.
    Each effect type has an associated impact level and can be mapped to a Sequence Ontology term.
    """
    # High impact effects
    CHROMOSOME_LARGE_DELETION = "CHROMOSOME_LARGE_DELETION"
    CHROMOSOME_LARGE_INVERSION = "CHROMOSOME_LARGE_INVERSION"
    CHROMOSOME_LARGE_DUPLICATION = "CHROMOSOME_LARGE_DUPLICATION"
    GENE_REARRANGEMENT = "GENE_REARRANGEMENT"
    GENE_DELETED = "GENE_DELETED"
    TRANSCRIPT_DELETED = "TRANSCRIPT_DELETED"
    EXON_DELETED = "EXON_DELETED"
    EXON_DELETED_PARTIAL = "EXON_DELETED_PARTIAL"
    GENE_FUSION = "GENE_FUSION"
    GENE_FUSION_REVERESE = "GENE_FUSION_REVERESE"
    GENE_FUSION_HALF = "GENE_FUSION_HALF"
    FRAME_SHIFT = "FRAME_SHIFT"
    STOP_GAINED = "STOP_GAINED"
    STOP_LOST = "STOP_LOST"
    START_LOST = "START_LOST"
    SPLICE_SITE_ACCEPTOR = "SPLICE_SITE_ACCEPTOR"
    SPLICE_SITE_DONOR = "SPLICE_SITE_DONOR"
    RARE_AMINO_ACID = "RARE_AMINO_ACID"
    EXON_DUPLICATION = "EXON_DUPLICATION"
    EXON_DUPLICATION_PARTIAL = "EXON_DUPLICATION_PARTIAL"
    EXON_INVERSION = "EXON_INVERSION"
    EXON_INVERSION_PARTIAL = "EXON_INVERSION_PARTIAL"
    PROTEIN_PROTEIN_INTERACTION_LOCUS = "PROTEIN_PROTEIN_INTERACTION_LOCUS"
    PROTEIN_STRUCTURAL_INTERACTION_LOCUS = "PROTEIN_STRUCTURAL_INTERACTION_LOCUS"
    
    # Moderate impact effects
    NON_SYNONYMOUS_CODING = "NON_SYNONYMOUS_CODING"
    GENE_DUPLICATION = "GENE_DUPLICATION"
    TRANSCRIPT_DUPLICATION = "TRANSCRIPT_DUPLICATION"
    UTR_5_DELETED = "UTR_5_DELETED"
    UTR_3_DELETED = "UTR_3_DELETED"
    SPLICE_SITE_BRANCH_U12 = "SPLICE_SITE_BRANCH_U12"
    GENE_INVERSION = "GENE_INVERSION"
    TRANSCRIPT_INVERSION = "TRANSCRIPT_INVERSION"
    CODON_INSERTION = "CODON_INSERTION"
    CODON_CHANGE_PLUS_CODON_INSERTION = "CODON_CHANGE_PLUS_CODON_INSERTION"
    CODON_DELETION = "CODON_DELETION"
    CODON_CHANGE_PLUS_CODON_DELETION = "CODON_CHANGE_PLUS_CODON_DELETION"
    
    # Low impact effects
    NON_SYNONYMOUS_STOP = "NON_SYNONYMOUS_STOP"
    NON_SYNONYMOUS_START = "NON_SYNONYMOUS_START"
    SPLICE_SITE_REGION = "SPLICE_SITE_REGION"
    SPLICE_SITE_BRANCH = "SPLICE_SITE_BRANCH"
    SYNONYMOUS_CODING = "SYNONYMOUS_CODING"
    SYNONYMOUS_START = "SYNONYMOUS_START"
    SYNONYMOUS_STOP = "SYNONYMOUS_STOP"
    CODON_CHANGE = "CODON_CHANGE"
    START_GAINED = "START_GAINED"
    MOTIF = "MOTIF"
    MOTIF_DELETED = "MOTIF_DELETED"
    FEATURE_FUSION = "FEATURE_FUSION"
    
    # Modifier impact effects
    FRAME_SHIFT_BEFORE_CDS_START = "FRAME_SHIFT_BEFORE_CDS_START"
    FRAME_SHIFT_AFTER_CDS_END = "FRAME_SHIFT_AFTER_CDS_END"
    UTR_5_PRIME = "UTR_5_PRIME"
    UTR_3_PRIME = "UTR_3_PRIME"
    REGULATION = "REGULATION"
    MICRO_RNA = "MICRO_RNA"
    UPSTREAM = "UPSTREAM"
    DOWNSTREAM = "DOWNSTREAM"
    NEXT_PROT = "NEXT_PROT"
    INTRON_CONSERVED = "INTRON_CONSERVED"
    INTRON = "INTRON"
    INTRAGENIC = "INTRAGENIC"
    INTERGENIC_CONSERVED = "INTERGENIC_CONSERVED"
    INTERGENIC = "INTERGENIC"
    CDS = "CDS"
    EXON = "EXON"
    TRANSCRIPT = "TRANSCRIPT"
    GENE = "GENE"
    SEQUENCE = "SEQUENCE"
    CHROMOSOME_ELONGATION = "CHROMOSOME_ELONGATION"
    CUSTOM = "CUSTOM"
    CHROMOSOME = "CHROMOSOME"
    GENOME = "GENOME"
    NONE = "NONE"
    
    def get_impact(self) -> AnnotationImpact:
        """
        Get the impact level for a given effect type.
        """
        if self in HIGH_IMPACT:
            return AnnotationImpact.HIGH
        elif self in MODERATE_IMPACT:
            return AnnotationImpact.MODERATE
        elif self in LOW_IMPACT:
            return AnnotationImpact.LOW
        else:
            return AnnotationImpact.MODIFIER
    
    def to_sequence_ontology(self, warn_on_missing: bool = False) -> str:
        """
        Convert the effect type to a Sequence Ontology term.
        """
        if self in EFFECT_TO_SEQUENCE_ONTOLOGY:
            return EFFECT_TO_SEQUENCE_ONTOLOGY[self]
        if warn_on_missing:
            import warnings
            warnings.warn(f"No Sequence Ontology mapping found for effect type: {self.value}")
        return ''

    @classmethod    
    def from_sequence_ontology(cls, so_term: str) -> "EffectType | None":
        """
        Convert a Sequence Ontology term to an EffectType.
        """        
        return SEQUENCE_ONTOLOGY_TO_EFFECT.get(so_term)
    

HIGH_IMPACT = {EffectType.CHROMOSOME_LARGE_DELETION, EffectType.CHROMOSOME_LARGE_INVERSION,
               EffectType.CHROMOSOME_LARGE_DUPLICATION, EffectType.GENE_REARRANGEMENT,
               EffectType.GENE_DELETED, EffectType.TRANSCRIPT_DELETED, EffectType.EXON_DELETED,
               EffectType.EXON_DELETED_PARTIAL, EffectType.GENE_FUSION, EffectType.GENE_FUSION_REVERESE,
               EffectType.GENE_FUSION_HALF, EffectType.FRAME_SHIFT, EffectType.STOP_GAINED,
               EffectType.STOP_LOST, EffectType.START_LOST, EffectType.SPLICE_SITE_ACCEPTOR,
               EffectType.SPLICE_SITE_DONOR, EffectType.RARE_AMINO_ACID, EffectType.EXON_DUPLICATION,
               EffectType.EXON_DUPLICATION_PARTIAL, EffectType.EXON_INVERSION, EffectType.EXON_INVERSION_PARTIAL,
               EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS, EffectType.PROTEIN_STRUCTURAL_INTERACTION_LOCUS
               }

MODERATE_IMPACT = { EffectType.NON_SYNONYMOUS_CODING, EffectType.GENE_DUPLICATION, EffectType.TRANSCRIPT_DUPLICATION,
                   EffectType.UTR_5_DELETED, EffectType.UTR_3_DELETED, EffectType.SPLICE_SITE_BRANCH_U12,
                   EffectType.GENE_INVERSION, EffectType.TRANSCRIPT_INVERSION, EffectType.CODON_INSERTION,
                   EffectType.CODON_CHANGE_PLUS_CODON_INSERTION, EffectType.CODON_DELETION,
                   EffectType.CODON_CHANGE_PLUS_CODON_DELETION
                   }
            
LOW_IMPACT = {EffectType.NON_SYNONYMOUS_STOP, EffectType.NON_SYNONYMOUS_START, EffectType.SPLICE_SITE_REGION,
              EffectType.SPLICE_SITE_BRANCH, EffectType.SYNONYMOUS_CODING, EffectType.SYNONYMOUS_START,
              EffectType.SYNONYMOUS_STOP, EffectType.CODON_CHANGE, EffectType.START_GAINED,
              EffectType.MOTIF, EffectType.MOTIF_DELETED, EffectType.FEATURE_FUSION
              }


EFFECT_TO_SEQUENCE_ONTOLOGY = {
    EffectType.CDS: "coding_sequence_variant",
    EffectType.CHROMOSOME_LARGE_DELETION: "chromosome_number_variation",
    EffectType.CHROMOSOME_LARGE_DUPLICATION: "duplication",
    EffectType.CHROMOSOME_LARGE_INVERSION: "inversion",
    EffectType.CHROMOSOME: "chromosome",
    EffectType.CHROMOSOME_ELONGATION: "feature_elongation",
    EffectType.CODON_CHANGE: "coding_sequence_variant",
    EffectType.CODON_CHANGE_PLUS_CODON_INSERTION: "disruptive_inframe_insertion",
    EffectType.CODON_CHANGE_PLUS_CODON_DELETION: "disruptive_inframe_deletion",
    EffectType.CODON_DELETION: "conservative_inframe_deletion",
    EffectType.CODON_INSERTION: "conservative_inframe_insertion",
    EffectType.DOWNSTREAM: "downstream_gene_variant",
    EffectType.EXON: "non_coding_transcript_exon_variant",
    EffectType.EXON_DELETED: "exon_loss_variant",
    EffectType.EXON_DELETED_PARTIAL: "exon_loss_variant",
    EffectType.EXON_DUPLICATION: "duplication",
    EffectType.EXON_DUPLICATION_PARTIAL: "duplication",
    EffectType.EXON_INVERSION: "inversion",
    EffectType.EXON_INVERSION_PARTIAL: "inversion",
    EffectType.FEATURE_FUSION: "feature_fusion",
    EffectType.FRAME_SHIFT: "frameshift_variant",
    EffectType.FRAME_SHIFT_BEFORE_CDS_START: "start_retained_variant",
    EffectType.FRAME_SHIFT_AFTER_CDS_END: "stop_retained_variant",
    EffectType.GENE: "gene_variant",
    EffectType.GENE_INVERSION: "inversion",
    EffectType.GENE_DELETED: "feature_ablation",
    EffectType.GENE_DUPLICATION: "duplication",
    EffectType.GENE_FUSION: "gene_fusion",
    EffectType.GENE_FUSION_HALF: "transcript_ablation",
    EffectType.GENE_FUSION_REVERESE: "bidirectional_gene_fusion",
    EffectType.GENE_REARRANGEMENT: "rearranged_at_DNA_level",
    EffectType.INTERGENIC: "intergenic_region",
    EffectType.INTERGENIC_CONSERVED: "conserved_intergenic_variant",
    EffectType.INTRON: "intron_variant",
    EffectType.INTRON_CONSERVED: "conserved_intron_variant",
    EffectType.INTRAGENIC: "intragenic_variant",
    EffectType.MICRO_RNA: "miRNA",
    EffectType.MOTIF: "TF_binding_site_variant",
    EffectType.MOTIF_DELETED: "TFBS_ablation",
    EffectType.NEXT_PROT: "sequence_feature",
    EffectType.NON_SYNONYMOUS_CODING: "missense_variant",
    EffectType.NON_SYNONYMOUS_START: "initiator_codon_variant",
    EffectType.NON_SYNONYMOUS_STOP: "stop_retained_variant",
    EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS: "protein_protein_contact",
    EffectType.PROTEIN_STRUCTURAL_INTERACTION_LOCUS: "structural_interaction_variant",
    EffectType.RARE_AMINO_ACID: "rare_amino_acid_variant",
    EffectType.REGULATION: "regulatory_region_variant",
    EffectType.SPLICE_SITE_ACCEPTOR: "splice_acceptor_variant",
    EffectType.SPLICE_SITE_DONOR: "splice_donor_variant",
    EffectType.SPLICE_SITE_REGION: "splice_region_variant",
    EffectType.SPLICE_SITE_BRANCH: "splice_branch_variant",
    EffectType.SPLICE_SITE_BRANCH_U12: "splice_branch_variant",
    EffectType.START_LOST: "start_lost",
    EffectType.START_GAINED: "5_prime_UTR_premature_start_codon_gain_variant",
    EffectType.STOP_GAINED: "stop_gained",
    EffectType.STOP_LOST: "stop_lost",
    EffectType.SYNONYMOUS_CODING: "synonymous_variant",
    EffectType.SYNONYMOUS_STOP: "stop_retained_variant",
    EffectType.SYNONYMOUS_START: "start_retained_variant",
    EffectType.TRANSCRIPT: "non_coding_transcript_variant",
    EffectType.TRANSCRIPT_DELETED: "transcript_ablation",
    EffectType.TRANSCRIPT_DUPLICATION: "duplication",
    EffectType.TRANSCRIPT_INVERSION: "inversion",
    EffectType.UPSTREAM: "upstream_gene_variant",
    EffectType.UTR_3_PRIME: "3_prime_UTR_variant",
    EffectType.UTR_3_DELETED: "3_prime_UTR_truncation+exon_loss_variant",
    EffectType.UTR_5_PRIME: "5_prime_UTR_variant",
    EffectType.UTR_5_DELETED: "5_prime_UTR_truncation+exon_loss_variant",
    EffectType.CUSTOM: "custom",
    }

SEQUENCE_ONTOLOGY_TO_EFFECT = {v: k for k, v in EFFECT_TO_SEQUENCE_ONTOLOGY.items()}