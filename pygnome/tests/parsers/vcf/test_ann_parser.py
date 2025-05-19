"""
Tests for the ANN field parser.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.ann import AnnParser, VariantAnnotation, AnnotationImpact, FeatureType
from pygnome.parsers.vcf.ann.effect_type import EffectType


class TestAnnParser(unittest.TestCase):
    """Test cases for the ANN field parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a VCF header with ANN field definition
        self.header = VcfHeader()
        self.header.add_meta_line('##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Putative_impact | Gene_name | Gene_ID | Feature_type | Feature_ID | Transcript_biotype | Rank | HGVS.c | HGVS.p | cDNA_position | CDS_position | Protein_position | Distance_to_feature | Errors_warnings_info\'">')
    
    def test_simple_annotation(self):
        """Test parsing a simple annotation."""
        # Create a VCF record with an ANN field
        record_line = "chr1\t100\t.\tA\tG\t.\t.\tANN=G|missense_variant|MODERATE|GENE1|ENSG00000123|transcript|ENST00000456|Coding|1/5|c.100A>G|p.Lys34Arg|100/1000|100/900|34/300||"
        record = VcfRecord(record_line, self.header)
        
        # Create a parser with the record
        parser = AnnParser(record)
        
        # Iterate over the annotations
        annotations = list(parser)
        
        # Check that we got one annotation
        self.assertEqual(len(annotations), 1)
        
        # Check the annotation fields
        ann = annotations[0]
        self.assertEqual(ann.allele, "G")
        self.assertEqual(ann.annotation, "missense_variant")
        self.assertEqual(ann.putative_impact, AnnotationImpact.MODERATE)
        self.assertEqual(ann.gene_name, "GENE1")
        self.assertEqual(ann.gene_id, "ENSG00000123")
        self.assertEqual(ann.feature_type, FeatureType.TRANSCRIPT)
        self.assertEqual(ann.feature_id, "ENST00000456")
        self.assertEqual(ann.rank, 1)
        self.assertEqual(ann.total, 5)
        self.assertEqual(ann.hgvs_c, "c.100A>G")
        self.assertEqual(ann.hgvs_p, "p.Lys34Arg")
        self.assertEqual(ann.cdna_pos, 100)
        self.assertEqual(ann.cdna_length, 1000)
        self.assertEqual(ann.cds_pos, 100)
        self.assertEqual(ann.cds_length, 900)
        self.assertEqual(ann.protein_pos, 34)
        self.assertEqual(ann.protein_length, 300)
        self.assertIsNone(ann.distance)
        self.assertIsNone(ann.messages)
    
    def test_multiple_annotations(self):
        """Test parsing multiple annotations."""
        # Create a VCF record with multiple annotations
        record_line = "chr1\t100\t.\tA\tG\t.\t.\tANN=G|missense_variant|MODERATE|GENE1|ENSG00000123|transcript|ENST00000456|Coding|1/5|c.100A>G|p.Lys34Arg|100/1000|100/900|34/300||,G|upstream_gene_variant|MODIFIER|GENE2|ENSG00000789|transcript|ENST00000789|Noncoding|||||||500|"
        record = VcfRecord(record_line, self.header)
        
        # Create a parser with the record
        parser = AnnParser(record)
        
        # Iterate over the annotations
        annotations = list(parser)
        
        # Check that we got two annotations
        self.assertEqual(len(annotations), 2)
        
        # Check the first annotation
        ann1 = annotations[0]
        self.assertEqual(ann1.allele, "G")
        self.assertEqual(ann1.annotation, "missense_variant")
        self.assertEqual(ann1.effect, EffectType.NON_SYNONYMOUS_CODING)
        self.assertEqual(ann1.putative_impact, AnnotationImpact.MODERATE)
        
        # Check the second annotation
        ann2 = annotations[1]
        self.assertEqual(ann2.allele, "G")
        self.assertEqual(ann2.annotation, "upstream_gene_variant")
        self.assertEqual(ann2.effect, EffectType.UPSTREAM)
        self.assertEqual(ann2.putative_impact, AnnotationImpact.MODIFIER)
        self.assertEqual(ann2.distance, 500)
    
    def test_with_errors(self):
        """Test parsing annotations with error codes."""
        # Create a VCF record with error codes
        record_line = "chr1\t100\t.\tA\tG\t.\t.\tANN=G|missense_variant|MODERATE|GENE1|ENSG00000123|transcript|ENST00000456|Coding|1/5|c.100A>G|p.Lys34Arg|100/1000|100/900|34/300||W1&W5"
        record = VcfRecord(record_line, self.header)
        
        # Create a parser with the record
        parser = AnnParser(record)
        
        # Iterate over the annotations
        annotations = list(parser)
        
        # Check that we got one annotation
        self.assertEqual(len(annotations), 1)
        
        # Check the error codes
        ann = annotations[0]
        self.assertIsNotNone(ann.messages)
        if ann.messages:  # Check if messages is not None
            self.assertEqual(len(ann.messages), 2)
            message_values = [m.value for m in ann.messages]
            self.assertIn("W1", message_values)
            self.assertIn("W5", message_values)
    
    def test_no_ann_field(self):
        """Test parsing a record without an ANN field."""
        # Create a VCF record without an ANN field
        record_line = "chr1\t100\t.\tA\tG\t.\t.\t."
        record = VcfRecord(record_line, self.header)
        
        # Create a parser with the record
        parser = AnnParser(record)
        
        # Iterate over the annotations
        annotations = list(parser)
        
        # Check that we got no annotations
        self.assertEqual(len(annotations), 0)
    
    def test_effect_type_mapping(self):
        """Test mapping between effect types and Sequence Ontology terms."""
        # Test mapping from effect type to SO term
        self.assertEqual(EffectType.NON_SYNONYMOUS_CODING.to_sequence_ontology(), "missense_variant")
        self.assertEqual(EffectType.FRAME_SHIFT.to_sequence_ontology(), "frameshift_variant")
        self.assertEqual(EffectType.STOP_GAINED.to_sequence_ontology(), "stop_gained")
        
        # Test mapping from SO term to effect type
        self.assertEqual(EffectType.from_sequence_ontology("missense_variant"), EffectType.NON_SYNONYMOUS_CODING)
        self.assertEqual(EffectType.from_sequence_ontology("frameshift_variant"), EffectType.FRAME_SHIFT)
        self.assertEqual(EffectType.from_sequence_ontology("stop_gained"), EffectType.STOP_GAINED)
    
    def test_effect_type_impact(self):
        """Test mapping between effect types and impact levels."""
        # Test high impact effects
        self.assertEqual(EffectType.FRAME_SHIFT.get_impact(), AnnotationImpact.HIGH)
        self.assertEqual(EffectType.STOP_GAINED.get_impact(), AnnotationImpact.HIGH)
        self.assertEqual(EffectType.SPLICE_SITE_ACCEPTOR.get_impact(), AnnotationImpact.HIGH)
        
        # Test moderate impact effects
        self.assertEqual(EffectType.NON_SYNONYMOUS_CODING.get_impact(), AnnotationImpact.MODERATE)
        self.assertEqual(EffectType.CODON_INSERTION.get_impact(), AnnotationImpact.MODERATE)
        
        # Test low impact effects
        self.assertEqual(EffectType.SYNONYMOUS_CODING.get_impact(), AnnotationImpact.LOW)
        self.assertEqual(EffectType.SPLICE_SITE_REGION.get_impact(), AnnotationImpact.LOW)
        
        # Test modifier impact effects
        self.assertEqual(EffectType.INTRON.get_impact(), AnnotationImpact.MODIFIER)
        self.assertEqual(EffectType.UPSTREAM.get_impact(), AnnotationImpact.MODIFIER)
    
    def test_effect_type_parsing(self):
        """Test parsing different annotation formats."""
        # Create a VCF record with an effect type instead of SO term
        record_line = "chr1\t100\t.\tA\tG\t.\t.\tANN=G|NON_SYNONYMOUS_CODING|MODERATE|GENE1|ENSG00000123|transcript|ENST00000456|Coding|1/5|c.100A>G|p.Lys34Arg|100/1000|100/900|34/300||"
        record = VcfRecord(record_line, self.header)
        
        # Create a parser with the record
        parser = AnnParser(record)
        
        # Iterate over the annotations
        annotations = list(parser)
        
        # Check that we got one annotation
        self.assertEqual(len(annotations), 1)
        
        # Check the annotation fields
        ann = annotations[0]
        self.assertEqual(ann.allele, "G")
        self.assertEqual(ann.annotation, "NON_SYNONYMOUS_CODING")
        self.assertEqual(ann.effect, EffectType.NON_SYNONYMOUS_CODING)
        self.assertEqual(ann.putative_impact, AnnotationImpact.MODERATE)


if __name__ == "__main__":
    unittest.main()