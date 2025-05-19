# To do

According to the 'VCF annotation format' document, the Effect types are defined ad have an associated impact level. 

Nevertheless in out `VariantAnnotation` class, we have `annotation: str`, which is clearly wrong.
Create an Enum with the effect types, and map them to `AnnotationImpact`.

Here is the Java code from SnpEff:

```java
public enum EffectType {
	// High impact
	// Order: Highest impact first
	CHROMOSOME_LARGE_DELETION(EffectImpact.HIGH) //
	, CHROMOSOME_LARGE_INVERSION(EffectImpact.HIGH) //
	, CHROMOSOME_LARGE_DUPLICATION(EffectImpact.HIGH) //
	, GENE_REARRANGEMENT(EffectImpact.HIGH) //
	, GENE_DELETED(EffectImpact.HIGH) //
	, TRANSCRIPT_DELETED(EffectImpact.HIGH) //
	, EXON_DELETED(EffectImpact.HIGH) //
	, EXON_DELETED_PARTIAL(EffectImpact.HIGH) //
	, GENE_FUSION(EffectImpact.HIGH) //
	, GENE_FUSION_REVERESE(EffectImpact.HIGH) //
	, GENE_FUSION_HALF(EffectImpact.HIGH) //
	, FRAME_SHIFT(EffectImpact.HIGH) //
	, STOP_GAINED(EffectImpact.HIGH) //
	, STOP_LOST(EffectImpact.HIGH) //
	, START_LOST(EffectImpact.HIGH) //
	, SPLICE_SITE_ACCEPTOR(EffectImpact.HIGH) //
	, SPLICE_SITE_DONOR(EffectImpact.HIGH) //
	, RARE_AMINO_ACID(EffectImpact.HIGH) //
	, EXON_DUPLICATION(EffectImpact.HIGH) //
	, EXON_DUPLICATION_PARTIAL(EffectImpact.HIGH) //
	, EXON_INVERSION(EffectImpact.HIGH) //
	, EXON_INVERSION_PARTIAL(EffectImpact.HIGH) //
	, PROTEIN_PROTEIN_INTERACTION_LOCUS(EffectImpact.HIGH) //
	, PROTEIN_STRUCTURAL_INTERACTION_LOCUS(EffectImpact.HIGH) //

	// Moderate impact
	// Order: Highest impact first
	// Note: Method Codon.effect() relies on this order for effect
	//       replacement (when 'allowReplace = true')
	, NON_SYNONYMOUS_CODING(EffectImpact.MODERATE) //
	, GENE_DUPLICATION(EffectImpact.MODERATE) //
	, TRANSCRIPT_DUPLICATION(EffectImpact.MODERATE) //
	, UTR_5_DELETED(EffectImpact.MODERATE) //
	, UTR_3_DELETED(EffectImpact.MODERATE) //
	, SPLICE_SITE_BRANCH_U12(EffectImpact.MODERATE) //
	, GENE_INVERSION(EffectImpact.MODERATE) //
	, TRANSCRIPT_INVERSION(EffectImpact.MODERATE) //
	, CODON_INSERTION(EffectImpact.MODERATE) //
	, CODON_CHANGE_PLUS_CODON_INSERTION(EffectImpact.MODERATE) //
	, CODON_DELETION(EffectImpact.MODERATE) //
	, CODON_CHANGE_PLUS_CODON_DELETION(EffectImpact.MODERATE) //

	// Low impact
	// Order: Highest impact first
	, NON_SYNONYMOUS_STOP(EffectImpact.LOW) //
	, NON_SYNONYMOUS_START(EffectImpact.LOW) //
	, SPLICE_SITE_REGION(EffectImpact.LOW) //
	, SPLICE_SITE_BRANCH(EffectImpact.LOW) //
	, SYNONYMOUS_CODING(EffectImpact.LOW) //
	, SYNONYMOUS_START(EffectImpact.LOW) //
	, SYNONYMOUS_STOP(EffectImpact.LOW) //
	, CODON_CHANGE(EffectImpact.LOW) //
	, START_GAINED(EffectImpact.LOW) //
	, MOTIF(EffectImpact.LOW) //
	, MOTIF_DELETED(EffectImpact.LOW) //
	, FEATURE_FUSION(EffectImpact.LOW) //

	// Modifiers
	// Order: Highest impact first
	, FRAME_SHIFT_BEFORE_CDS_START(EffectImpact.MODIFIER) //
	, FRAME_SHIFT_AFTER_CDS_END(EffectImpact.MODIFIER) //
	, UTR_5_PRIME(EffectImpact.MODIFIER) //
	, UTR_3_PRIME(EffectImpact.MODIFIER) //
	, REGULATION(EffectImpact.MODIFIER) //
	, MICRO_RNA(EffectImpact.MODIFIER) //
	, UPSTREAM(EffectImpact.MODIFIER) //
	, DOWNSTREAM(EffectImpact.MODIFIER) //
	, NEXT_PROT(EffectImpact.MODIFIER) //
	, INTRON_CONSERVED(EffectImpact.MODIFIER) //
	, INTRON(EffectImpact.MODIFIER) //
	, INTRAGENIC(EffectImpact.MODIFIER) //
	, INTERGENIC_CONSERVED(EffectImpact.MODIFIER) //
	, INTERGENIC(EffectImpact.MODIFIER) //
	, CDS(EffectImpact.MODIFIER) //
	, EXON(EffectImpact.MODIFIER) //
	, TRANSCRIPT(EffectImpact.MODIFIER) //
	, GENE(EffectImpact.MODIFIER) //
	, SEQUENCE(EffectImpact.MODIFIER) //
	, CHROMOSOME_ELONGATION(EffectImpact.MODIFIER) //
	, CUSTOM(EffectImpact.MODIFIER) //
	, CHROMOSOME(EffectImpact.MODIFIER) //
	, GENOME(EffectImpact.MODIFIER) //
	, NONE(EffectImpact.MODIFIER) //
	;
    ...
}
```

Also these Effect types should be represented as "Sequence Ontology" terms.
Here is the mapping (also from SnpEff):

```java
	public String toSequenceOntology(EffFormatVersion formatVersion, Variant variant) {
		switch (this) {

		case CDS:
			return "coding_sequence_variant";

		case CHROMOSOME_LARGE_DELETION:
			return "chromosome_number_variation";

		case CHROMOSOME_LARGE_DUPLICATION:
			return "duplication";

		case CHROMOSOME_LARGE_INVERSION:
			return "inversion";

		case CHROMOSOME:
			return "chromosome";

		case CHROMOSOME_ELONGATION:
			return "feature_elongation";

		case CODON_CHANGE:
			return "coding_sequence_variant";

		case CODON_CHANGE_PLUS_CODON_INSERTION:
			return "disruptive_inframe_insertion";

		case CODON_CHANGE_PLUS_CODON_DELETION:
			return "disruptive_inframe_deletion";

		case CODON_DELETION:
			return "conservative_inframe_deletion";

		case CODON_INSERTION:
			return "conservative_inframe_insertion";

		case DOWNSTREAM:
			return "downstream_gene_variant";

		case EXON:
			if (variant != null && (!variant.isVariant() || variant.isInterval())) return "exon_region";
			return "non_coding_transcript_exon_variant";

		case EXON_DELETED:
			return "exon_loss_variant";

		case EXON_DELETED_PARTIAL:
			return "exon_loss_variant";

		case EXON_DUPLICATION:
			return "duplication";

		case EXON_DUPLICATION_PARTIAL:
			return "duplication";

		case EXON_INVERSION:
			return "inversion";

		case EXON_INVERSION_PARTIAL:
			return "inversion";

		case FEATURE_FUSION:
			return "feature_fusion";

		case FRAME_SHIFT:
			return "frameshift_variant";

		case FRAME_SHIFT_BEFORE_CDS_START:
			return "start_retained_variant";

		case FRAME_SHIFT_AFTER_CDS_END:
			return "stop_retained_variant";

		case GENE:
			return "gene_variant";

		case GENE_INVERSION:
			return "inversion";

		case GENE_DELETED:
			return "feature_ablation";

		case GENE_DUPLICATION:
			return "duplication";

		case GENE_FUSION:
			return "gene_fusion";

		case GENE_FUSION_HALF:
			return "transcript_ablation";

		case GENE_FUSION_REVERESE:
			return "bidirectional_gene_fusion";

		case GENE_REARRANGEMENT:
			return "rearranged_at_DNA_level";

		case INTERGENIC:
			return "intergenic_region";

		case INTERGENIC_CONSERVED:
			return "conserved_intergenic_variant";

		case INTRON:
			return "intron_variant";

		case INTRON_CONSERVED:
			return "conserved_intron_variant";

		case INTRAGENIC:
			return "intragenic_variant";

		case MICRO_RNA:
			return "miRNA";

		case MOTIF:
			return "TF_binding_site_variant";

		case MOTIF_DELETED:
			return "TFBS_ablation";

		case NEXT_PROT:
			return "sequence_feature";

		case NON_SYNONYMOUS_CODING:
			return "missense_variant";

		case NON_SYNONYMOUS_START:
			return "initiator_codon_variant";

		case NON_SYNONYMOUS_STOP:
			return "stop_retained_variant";

		case PROTEIN_PROTEIN_INTERACTION_LOCUS:
			return "protein_protein_contact";

		case PROTEIN_STRUCTURAL_INTERACTION_LOCUS:
			return "structural_interaction_variant";

		case RARE_AMINO_ACID:
			return "rare_amino_acid_variant";

		case REGULATION:
			return "regulatory_region_variant";

		case SPLICE_SITE_ACCEPTOR:
			return "splice_acceptor_variant";

		case SPLICE_SITE_DONOR:
			return "splice_donor_variant";

		case SPLICE_SITE_REGION:
			return "splice_region_variant";

		case SPLICE_SITE_BRANCH:
			return "splice_branch_variant";

		case SPLICE_SITE_BRANCH_U12:
			return "splice_branch_variant";

		case START_LOST:
			return "start_lost";

		case START_GAINED:
			return "5_prime_UTR_premature_start_codon_gain_variant";

		case STOP_GAINED:
			return "stop_gained";

		case STOP_LOST:
			return "stop_lost";

		case SYNONYMOUS_CODING:
			return "synonymous_variant";

		case SYNONYMOUS_STOP:
			return "stop_retained_variant";

		case SYNONYMOUS_START:
			return "start_retained_variant";

		case TRANSCRIPT:
			return "non_coding_transcript_variant";

		case TRANSCRIPT_DELETED:
			return "transcript_ablation";

		case TRANSCRIPT_DUPLICATION:
			return "duplication";

		case TRANSCRIPT_INVERSION:
			return "inversion";

		case UPSTREAM:
			return "upstream_gene_variant";

		case UTR_3_PRIME:
			return "3_prime_UTR_variant";

		case UTR_3_DELETED:
			return "3_prime_UTR_truncation" + formatVersion.separator() + "exon_loss_variant";

		case UTR_5_PRIME:
			return "5_prime_UTR_variant";

		case UTR_5_DELETED:
			return "5_prime_UTR_truncation" + formatVersion.separator() + "exon_loss_variant";

		case CUSTOM:
			return "custom";

		case NONE:
		case GENOME:
		case SEQUENCE:
			return "";

		default:
			throw new RuntimeException("Unknown SO term for " + this.toString());
		}
	}
```


1) Create the enums and mapping methods.
2) When parsing the VCF, typically you'd receive "Sequence ontology" terms, but you may receive 'Effect Type' terms as well. 
3) Internally we'll use EffectTypes, but when writing the VCF, we should use the Sequence Ontology terms.
