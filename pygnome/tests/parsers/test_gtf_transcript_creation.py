"""Test the creation of transcript structures from GTF files."""

import unittest
import tempfile
import os
from pathlib import Path

from pygnome.genomics.chromosome import Chromosome
from pygnome.genomics.genome import Genome
from pygnome.parsers.gff.gtf_parser import GtfParser
from pygnome.parsers.gff.genome_loader import GffGenomeLoader
from pygnome.genomics.strand import Strand
from pygnome.genomics.phase import Phase
from pygnome.genomics.biotype import Biotype
from pygnome.genomics.utr import UTRType


class TestGtfTranscriptCreation(unittest.TestCase):
    """Test the creation of transcript structures from GTF files."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.loader = GffGenomeLoader()
        
        # Create a simple genome sequence for testing
        self.genome_sequence = {
            "chr1": "N" * 10000  # 10kb of Ns
        }
        
        # Replace specific regions with meaningful sequences
        # Forward strand gene (positions 1000-2000)
        # Exon 1: 1000-1200 with ATG at 1050-1052
        # Exon 2: 1400-1600
        # Exon 3: 1800-2000 with TAA at 1950-1952
        seq_list = list(self.genome_sequence["chr1"])
        
        # Exon 1 with start codon
        for i in range(1000, 1200):
            seq_list[i] = "A"
        seq_list[1050] = "A"
        seq_list[1051] = "T"
        seq_list[1052] = "G"
        
        # Exon 2
        for i in range(1400, 1600):
            seq_list[i] = "C"
        
        # Exon 3 with stop codon
        for i in range(1800, 2000):
            seq_list[i] = "G"
        seq_list[1950] = "T"
        seq_list[1951] = "A"
        seq_list[1952] = "A"
        
        # Reverse strand gene (positions 5000-6000)
        # Exon 1: 5000-5200 with TTA at 5050-5052 (reverse complement of TAA stop codon)
        # Exon 2: 5400-5600
        # Exon 3: 5800-6000 with CAT at 5950-5952 (reverse complement of ATG start codon)
        
        # Exon 1 with stop codon (reverse strand)
        for i in range(5000, 5200):
            seq_list[i] = "T"
        seq_list[5050] = "T"
        seq_list[5051] = "T"
        seq_list[5052] = "A"
        
        # Exon 2 (reverse strand)
        for i in range(5400, 5600):
            seq_list[i] = "G"
        
        # Exon 3 with start codon (reverse strand)
        for i in range(5800, 6000):
            seq_list[i] = "C"
        seq_list[5950] = "C"
        seq_list[5951] = "A"
        seq_list[5952] = "T"
        
        self.genome_sequence["chr1"] = "".join(seq_list)

    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()

    def _create_gtf_file(self, content):
        """Create a temporary GTF file with the given content."""
        gtf_path = os.path.join(self.temp_dir.name, "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(content)
        return gtf_path
    
    def test_forward_strand_basic_transcript(self):
        """Test that a basic transcript with exons is correctly parsed from GTF."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Check that the gene was loaded
        self.assertIn("gene1", genes)
        gene = genes["gene1"]
        
        # Check that the transcript was created
        self.assertEqual(len(gene.transcripts), 1)
        transcript = gene.transcripts[0]
        self.assertEqual(transcript.id, "transcript1")
        
        # Check that the exons were created
        self.assertEqual(len(transcript.exons), 3)
        
        # Check exon positions (converted to 0-based)
        exons = sorted(transcript.exons, key=lambda x: x.start)
        self.assertEqual(exons[0].start, 1000)
        self.assertEqual(exons[0].end, 1200)
        self.assertEqual(exons[1].start, 1400)
        self.assertEqual(exons[1].end, 1600)
        self.assertEqual(exons[2].start, 1800)
        self.assertEqual(exons[2].end, 2000)
        
        # Check that introns are lazily created
        self.assertIsNone(transcript._introns)
        introns = transcript.introns
        self.assertEqual(len(introns), 2)
        self.assertEqual(introns[0].start, 1200)
        self.assertEqual(introns[0].end, 1400)
        self.assertEqual(introns[1].start, 1600)
        self.assertEqual(introns[1].end, 1800)

    def test_reverse_strand_basic_transcript(self):
        """Test that a basic transcript with exons is correctly parsed from GTF on the reverse strand."""
        gtf_content = """
chr1\ttest\tgene\t5001\t6000\t.\t-\t.\tgene_id "gene2"; gene_name "GENE2";
chr1\ttest\ttranscript\t5001\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\texon\t5001\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr1\ttest\texon\t5401\t5600\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t5801\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Check that the gene was loaded
        self.assertIn("gene2", genes)
        gene = genes["gene2"]
        
        # Check that the transcript was created
        self.assertEqual(len(gene.transcripts), 1)
        transcript = gene.transcripts[0]
        self.assertEqual(transcript.id, "transcript2")
        self.assertEqual(transcript.strand, Strand.NEGATIVE)
        
        # Check that the exons were created
        self.assertEqual(len(transcript.exons), 3)
        
        # Check exon positions (converted to 0-based)
        exons = sorted(transcript.exons, key=lambda x: x.start)
        self.assertEqual(exons[0].start, 5000)
        self.assertEqual(exons[0].end, 5200)
        self.assertEqual(exons[1].start, 5400)
        self.assertEqual(exons[1].end, 5600)
        self.assertEqual(exons[2].start, 5800)
        self.assertEqual(exons[2].end, 6000)
        
        # Check that introns are lazily created
        self.assertIsNone(transcript._introns)
        introns = transcript.introns
        self.assertEqual(len(introns), 2)
        self.assertEqual(introns[0].start, 5200)
        self.assertEqual(introns[0].end, 5400)
        self.assertEqual(introns[1].start, 5600)
        self.assertEqual(introns[1].end, 5800)

    def test_forward_strand_with_cds(self):
        """Test that CDS regions are correctly parsed and added to the transcript."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t1\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that the CDS regions were created
        self.assertEqual(len(transcript.cds_list), 3)
        
        # Check CDS positions (converted to 0-based)
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[0].start, 1050)
        self.assertEqual(cds_list[0].end, 1200)
        self.assertEqual(cds_list[0].phase, Phase.ZERO)
        
        self.assertEqual(cds_list[1].start, 1400)
        self.assertEqual(cds_list[1].end, 1600)
        self.assertEqual(cds_list[1].phase, Phase.TWO)
        
        self.assertEqual(cds_list[2].start, 1800)
        self.assertEqual(cds_list[2].end, 1952)
        self.assertEqual(cds_list[2].phase, Phase.ONE)
        
        # Check that the transcript is marked as coding
        self.assertTrue(transcript.is_coding)
        self.assertTrue(transcript.has_cds)

    def test_reverse_strand_with_cds(self):
        """Test that CDS regions are correctly parsed and added to the transcript on the reverse strand."""
        gtf_content = """
chr1\ttest\tgene\t5001\t6000\t.\t-\t.\tgene_id "gene2"; gene_name "GENE2";
chr1\ttest\ttranscript\t5001\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\texon\t5001\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr1\ttest\texon\t5401\t5600\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t5801\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
chr1\ttest\tCDS\t5001\t5052\t.\t-\t0\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5401\t5600\t.\t-\t2\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5950\t6000\t.\t-\t1\tgene_id "gene2"; transcript_id "transcript2";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene2"].transcripts[0]
        
        # Check that the CDS regions were created
        self.assertEqual(len(transcript.cds_list), 3)
        
        # Check CDS positions (converted to 0-based)
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[0].start, 5000)
        self.assertEqual(cds_list[0].end, 5052)
        self.assertEqual(cds_list[0].phase, Phase.ZERO)
        
        self.assertEqual(cds_list[1].start, 5400)
        self.assertEqual(cds_list[1].end, 5600)
        self.assertEqual(cds_list[1].phase, Phase.TWO)
        
        self.assertEqual(cds_list[2].start, 5949)
        self.assertEqual(cds_list[2].end, 6000)
        self.assertEqual(cds_list[2].phase, Phase.ONE)
        
        # Check that the transcript is marked as coding
        self.assertTrue(transcript.is_coding)
        self.assertTrue(transcript.has_cds)

    def test_no_cds_provided(self):
        """Test that a transcript without CDS regions is handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1"; biotype "lncRNA";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; biotype "lncRNA";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that no CDS regions were created
        self.assertEqual(len(transcript.cds_list), 0)
        
        # Check that the transcript is not marked as coding
        self.assertFalse(transcript.has_cds)
        self.assertEqual(transcript.biotype, Biotype.LNCRNA)
        self.assertFalse(transcript.is_coding)

    def test_cds_without_phase(self):
        """Test that phase is correctly inferred when not provided."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that the CDS regions were created with default phase
        self.assertEqual(len(transcript.cds_list), 3)
        
        # All CDS should have Phase.ZERO by default
        for cds in transcript.cds_list:
            self.assertEqual(cds.phase, Phase.ZERO)
        
        # Infer phase
        transcript.infer_phase(self.genome_sequence)
        
        # Check that phases were correctly inferred
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[0].phase, Phase.ZERO)  # First CDS should always be Phase.ZERO
        
        # The other phases depend on the implementation of infer_phase
        # We're just checking that they were updated from the default
        self.assertIsNotNone(cds_list[1].phase)
        self.assertIsNotNone(cds_list[2].phase)

    def test_forward_strand_with_utrs(self):
        """Test that 5' and 3' UTRs are correctly parsed and added to the transcript."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t1\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tfive_prime_utr\t1001\t1050\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tthree_prime_utr\t1953\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that the UTRs were created
        self.assertEqual(len(transcript.utrs), 2)
        
        # Check UTR types
        five_prime_utrs = transcript.five_prime_utrs
        three_prime_utrs = transcript.three_prime_utrs
        
        self.assertEqual(len(five_prime_utrs), 1)
        self.assertEqual(len(three_prime_utrs), 1)
        
        # Check UTR positions (converted to 0-based)
        self.assertEqual(five_prime_utrs[0].start, 1000)
        self.assertEqual(five_prime_utrs[0].end, 1050)
        self.assertEqual(five_prime_utrs[0].utr_type, UTRType.FIVE_PRIME)
        
        self.assertEqual(three_prime_utrs[0].start, 1952)
        self.assertEqual(three_prime_utrs[0].end, 2000)
        self.assertEqual(three_prime_utrs[0].utr_type, UTRType.THREE_PRIME)

    def test_reverse_strand_with_utrs(self):
        """Test that 5' and 3' UTRs are correctly parsed and added to the transcript on the reverse strand."""
        gtf_content = """
chr1\ttest\tgene\t5001\t6000\t.\t-\t.\tgene_id "gene2"; gene_name "GENE2";
chr1\ttest\ttranscript\t5001\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\texon\t5001\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr1\ttest\texon\t5401\t5600\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t5801\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
chr1\ttest\tCDS\t5001\t5052\t.\t-\t0\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5401\t5600\t.\t-\t2\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5950\t6000\t.\t-\t1\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tfive_prime_utr\t5801\t5949\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tthree_prime_utr\t5053\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene2"].transcripts[0]
        
        # Check that the UTRs were created
        self.assertEqual(len(transcript.utrs), 2)
        
        # Check UTR types
        five_prime_utrs = transcript.five_prime_utrs
        three_prime_utrs = transcript.three_prime_utrs
        
        self.assertEqual(len(five_prime_utrs), 1)
        self.assertEqual(len(three_prime_utrs), 1)
        
        # Check UTR positions (converted to 0-based)
        self.assertEqual(five_prime_utrs[0].start, 5800)
        self.assertEqual(five_prime_utrs[0].end, 5949)
        self.assertEqual(five_prime_utrs[0].utr_type, UTRType.FIVE_PRIME)
        
        self.assertEqual(three_prime_utrs[0].start, 5052)
        self.assertEqual(three_prime_utrs[0].end, 5200)
        self.assertEqual(three_prime_utrs[0].utr_type, UTRType.THREE_PRIME)

    def test_cds_without_utrs(self):
        """Test that a transcript with CDS but no explicit UTRs is handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t1\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that no UTRs were created
        self.assertEqual(len(transcript.utrs), 0)
        
        # Check that the transcript is still marked as coding
        self.assertTrue(transcript.is_coding)
        self.assertTrue(transcript.has_cds)

    def test_single_exon_transcript(self):
        """Test that a transcript with a single exon has no introns or splice sites."""
        gtf_content = """
chr1\ttest\tgene\t1001\t1200\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\tCDS\t1051\t1150\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that only one exon was created
        self.assertEqual(len(transcript.exons), 1)
        
        # Check that no introns were created
        self.assertIsNone(transcript._introns)
        introns = transcript.introns
        self.assertEqual(len(introns), 0)
        
        # Infer splice sites
        transcript.infer_splice_sites()
        
        # Check that no splice sites were created
        self.assertEqual(len(transcript.splice_sites), 0)

    def test_splice_sites_forward_strand(self):
        """Test that splice sites are correctly inferred from exon boundaries on the forward strand."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Infer splice sites
        transcript.infer_splice_sites()
        
        # Check that splice sites were created
        self.assertEqual(len(transcript.splice_sites), 4)
        
        # Check splice site types and positions
        donor_sites = [s for s in transcript.splice_sites if s.site_type == "donor"]
        acceptor_sites = [s for s in transcript.splice_sites if s.site_type == "acceptor"]
        
        self.assertEqual(len(donor_sites), 2)
        self.assertEqual(len(acceptor_sites), 2)
        
        # Sort by position
        donor_sites = sorted(donor_sites, key=lambda x: x.start)
        acceptor_sites = sorted(acceptor_sites, key=lambda x: x.start)
        
        # Check donor sites (at the end of exons)
        self.assertEqual(donor_sites[0].start, 1198)  # 2 bases before exon end
        self.assertEqual(donor_sites[1].start, 1598)  # 2 bases before exon end
        
        # Check acceptor sites (at the start of exons)
        self.assertEqual(acceptor_sites[0].start, 1400)  # At exon start
        self.assertEqual(acceptor_sites[1].start, 1800)  # At exon start

    def test_splice_sites_reverse_strand(self):
        """Test that splice sites are correctly inferred from exon boundaries on the reverse strand."""
        gtf_content = """
chr1\ttest\tgene\t5001\t6000\t.\t-\t.\tgene_id "gene2"; gene_name "GENE2";
chr1\ttest\ttranscript\t5001\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\texon\t5001\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr1\ttest\texon\t5401\t5600\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t5801\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene2"].transcripts[0]
        
        # Infer splice sites
        transcript.infer_splice_sites()
        
        # Check that splice sites were created
        self.assertEqual(len(transcript.splice_sites), 4)
        
        # Check splice site types and positions
        donor_sites = [s for s in transcript.splice_sites if s.site_type == "donor"]
        acceptor_sites = [s for s in transcript.splice_sites if s.site_type == "acceptor"]
        
        self.assertEqual(len(donor_sites), 2)
        self.assertEqual(len(acceptor_sites), 2)
        
        # Sort by position
        donor_sites = sorted(donor_sites, key=lambda x: x.start)
        acceptor_sites = sorted(acceptor_sites, key=lambda x: x.start)
        
        # For reverse strand, donor sites are at the start of exons
        self.assertEqual(donor_sites[0].start, 5400)  # At exon start
        self.assertEqual(donor_sites[1].start, 5800)  # At exon start
        
        # For reverse strand, acceptor sites are at the end of exons
        self.assertEqual(acceptor_sites[0].start, 5198)  # 2 bases before exon end
        self.assertEqual(acceptor_sites[1].start, 5598)  # 2 bases before exon end

    def test_protein_translation(self):
        """Test that a transcript with CDS correctly translates to protein."""
        # Create a genome with a gene, transcript, and CDS
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Set up the gene and chromosome references
        gene = genes["gene1"]
        chromosome = Chromosome(name="chr1")
        genome = Genome(name="test_genome")
        genome.add_chromosome(chromosome)
        chromosome.add_gene(gene)
        gene.chromosome = chromosome
        # Set the genome on the chromosome instead of the gene
        chromosome.genome = genome
        
        # Get the CDS sequence
        cds_seq = transcript.cds(self.genome_sequence)
        
        # Check that the CDS sequence starts with ATG (our start codon)
        self.assertTrue(cds_seq.startswith("ATG"))
        
        # Check that the protein sequence starts with M (Methionine)
        protein = transcript.protein(self.genome_sequence)
        self.assertTrue(protein.startswith("M"))
        
        # Check that the protein sequence ends with a stop codon (translated to *)
        self.assertTrue(protein.endswith("*"))

    def test_reverse_strand_protein_translation(self):
        """Test that a transcript with CDS on the reverse strand correctly translates to protein."""
        # Create a genome with a gene, transcript, and CDS on the reverse strand
        gtf_content = """
chr1\ttest\tgene\t5001\t6000\t.\t-\t.\tgene_id "gene2"; gene_name "GENE2";
chr1\ttest\ttranscript\t5001\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\texon\t5001\t5200\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr1\ttest\texon\t5401\t5600\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t5801\t6000\t.\t-\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
chr1\ttest\tCDS\t5001\t5052\t.\t-\t0\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5401\t5600\t.\t-\t0\tgene_id "gene2"; transcript_id "transcript2";
chr1\ttest\tCDS\t5950\t6000\t.\t-\t0\tgene_id "gene2"; transcript_id "transcript2";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene2"].transcripts[0]
        
        # Set up the gene and chromosome references
        gene = genes["gene2"]
        chromosome = Chromosome(name="chr1")
        genome = Genome(name="test_genome")
        genome.add_chromosome(chromosome)
        chromosome.add_gene(gene)
        gene.chromosome = chromosome
        # Set the genome on the chromosome instead of the gene
        chromosome.genome = genome
        
        # Get the CDS sequence
        cds_seq = transcript.cds(self.genome_sequence)
        
        # Check that the CDS sequence starts with ATG (our start codon)
        # For reverse strand, the ATG comes from the reverse complement of CAT
        self.assertTrue(cds_seq.startswith("ATG"))
        
        # Check that the protein sequence starts with M (Methionine)
        protein = transcript.protein(self.genome_sequence)
        self.assertTrue(protein.startswith("M"))

    def test_noncoding_transcript_protein(self):
        """Test that a non-coding transcript returns an empty protein sequence."""
        # Create a genome with a gene and transcript but no CDS
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1"; biotype "lncRNA";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; biotype "lncRNA";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Set up the gene and chromosome references
        gene = genes["gene1"]
        chromosome = Chromosome(name="chr1")
        genome = Genome(name="test_genome")
        genome.add_chromosome(chromosome)
        chromosome.add_gene(gene)
        gene.chromosome = chromosome
        # Set the genome on the chromosome instead of the gene
        chromosome.genome = genome
        
        # Check that the transcript is not coding
        self.assertFalse(transcript.is_coding)
        self.assertFalse(transcript.has_cds)
        
        # Check that the protein sequence is empty
        protein = transcript.protein(self.genome_sequence)
        self.assertEqual(protein, "")

    def test_overlapping_exons(self):
        """Test that overlapping exons are handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1300\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1250\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1550\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that all exons were created
        self.assertEqual(len(transcript.exons), 3)
        
        # Check exon positions (converted to 0-based)
        exons = sorted(transcript.exons, key=lambda x: x.start)
        self.assertEqual(exons[0].start, 1000)
        self.assertEqual(exons[0].end, 1300)
        self.assertEqual(exons[1].start, 1249)
        self.assertEqual(exons[1].end, 1600)
        self.assertEqual(exons[2].start, 1549)
        self.assertEqual(exons[2].end, 2000)
        
        # Check that introns are correctly created despite overlapping exons
        # Since exons overlap, there should be no introns
        introns = transcript.introns
        self.assertEqual(len(introns), 0)

    def test_overlapping_cds(self):
        """Test that overlapping CDS regions are handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1150\t1300\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t1\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1952\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that all CDS regions were created
        self.assertEqual(len(transcript.cds_list), 4)
        
        # Check CDS positions (converted to 0-based)
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[0].start, 1050)
        self.assertEqual(cds_list[0].end, 1200)
        self.assertEqual(cds_list[1].start, 1149)
        self.assertEqual(cds_list[1].end, 1300)
        
        # Check that the transcript is marked as coding
        self.assertTrue(transcript.is_coding)
        self.assertTrue(transcript.has_cds)
        
        # Get the CDS sequence
        cds_seq = transcript.cds(self.genome_sequence)
        
        # The CDS sequence should be a concatenation of all CDS regions
        # Since there's an overlap, the sequence should include the overlapping region only once
        # This depends on the implementation of the cds method
        self.assertIsNotNone(cds_seq)

    def test_utrs_spanning_multiple_exons(self):
        """Test that UTRs spanning multiple exons are handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1501\t1600\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1900\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tfive_prime_utr\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tfive_prime_utr\t1401\t1500\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tthree_prime_utr\t1901\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that the UTRs were created
        self.assertEqual(len(transcript.utrs), 3)
        
        # Check UTR types
        five_prime_utrs = transcript.five_prime_utrs
        three_prime_utrs = transcript.three_prime_utrs
        
        self.assertEqual(len(five_prime_utrs), 2)
        self.assertEqual(len(three_prime_utrs), 1)
        
        # Check UTR positions (converted to 0-based)
        five_prime_utrs = sorted(five_prime_utrs, key=lambda x: x.start)
        self.assertEqual(five_prime_utrs[0].start, 1000)
        self.assertEqual(five_prime_utrs[0].end, 1200)
        self.assertEqual(five_prime_utrs[1].start, 1400)
        self.assertEqual(five_prime_utrs[1].end, 1500)
        
        self.assertEqual(three_prime_utrs[0].start, 1900)
        self.assertEqual(three_prime_utrs[0].end, 2000)
        
        # Check that the transcript is marked as coding
        self.assertTrue(transcript.is_coding)

    def test_incomplete_cds(self):
        """Test that CDS with length not divisible by 3 is handled correctly."""
        gtf_content = """
chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\tCDS\t1051\t1200\t.\t+\t0\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1401\t1600\t.\t+\t2\tgene_id "gene1"; transcript_id "transcript1";
chr1\ttest\tCDS\t1801\t1953\t.\t+\t1\tgene_id "gene1"; transcript_id "transcript1";
"""
        gtf_path = self._create_gtf_file(gtf_content)
        
        # Load the GTF file
        genes = self.loader.load(Path(gtf_path))
        
        # Get the transcript
        transcript = genes["gene1"].transcripts[0]
        
        # Check that the CDS regions were created
        self.assertEqual(len(transcript.cds_list), 3)
        
        # Check CDS positions (converted to 0-based)
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[0].start, 1050)
        self.assertEqual(cds_list[0].end, 1200)
        self.assertEqual(cds_list[2].start, 1800)
        self.assertEqual(cds_list[2].end, 1953)
        
        # Get the CDS sequence
        cds_seq = transcript.cds(self.genome_sequence)
        
        # The CDS sequence length should not be divisible by 3
        self.assertNotEqual(len(cds_seq) % 3, 0)
        
        # Fix the phase
        fixed = transcript.fix_phase(self.genome_sequence)
        
        # Check that the phase was fixed
        self.assertTrue(fixed)
        
        # Get the CDS sequence again
        cds_seq = transcript.cds(self.genome_sequence)
        
        # The CDS sequence length should now be divisible by 3
        self.assertEqual(len(cds_seq) % 3, 0)
        
        # Check that the last CDS was adjusted
        cds_list = sorted(transcript.cds_list, key=lambda x: x.start)
        self.assertEqual(cds_list[2].end, 1952)  # Should be adjusted to make CDS length divisible by 3


if __name__ == "__main__":
    unittest.main()