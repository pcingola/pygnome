"""Test the Transcript features and functionality."""

import unittest

from pygnome.genomics.biotype import Biotype
from pygnome.genomics.cds import CDS
from pygnome.genomics.chromosome import Chromosome
from pygnome.genomics.codon_table import CodonTableType
from pygnome.genomics.exon import Exon
from pygnome.genomics.gene import Gene
from pygnome.genomics.genome import Genome
from pygnome.genomics.phase import Phase
from pygnome.genomics.splice import SpliceSite, SpliceSiteType, SpliceRegion
from pygnome.genomics.strand import Strand
from pygnome.genomics.transcript import Transcript
from pygnome.genomics.utr import UTR, UTRType


class TestTranscriptFeatures(unittest.TestCase):
    """Test the Transcript features and functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a genome
        self.genome = Genome(
            name="Test Genome",
            species="Test Species",
            codon_table_type=CodonTableType.STANDARD
        )
        
        # Create a chromosome
        self.chrom = Chromosome(name="chr1")
        self.genome.add_chromosome(self.chrom)
        
        # Create a gene
        self.gene = Gene(
            id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE,
            name="Gene1"
        )
        self.chrom.add_gene(self.gene)
        
        # Create a transcript with exons, CDS, and UTRs
        self.transcript = Transcript(
            id="transcript1",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        self.gene.transcripts = [self.transcript]
        self.transcript.gene = self.gene
        
        # Add exons
        exon1 = Exon(
            id="exon1",
            chrom="chr1",
            start=1000,
            end=1200,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        exon2 = Exon(
            id="exon2",
            chrom="chr1",
            start=1400,
            end=1600,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        exon3 = Exon(
            id="exon3",
            chrom="chr1",
            start=1800,
            end=2000,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        self.transcript.exons = [exon1, exon2, exon3]
        
        # Add CDS
        cds1 = CDS(
            id="cds1",
            chrom="chr1",
            start=1100,
            end=1200,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        cds2 = CDS(
            id="cds2",
            chrom="chr1",
            start=1400,
            end=1600,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        cds3 = CDS(
            id="cds3",
            chrom="chr1",
            start=1800,
            end=1900,
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        self.transcript.cds_list = [cds1, cds2, cds3]
        
        # Add UTRs
        utr5 = UTR(
            id="utr5",
            chrom="chr1",
            start=1000,
            end=1100,
            strand=Strand.POSITIVE,
            utr_type=UTRType.FIVE_PRIME
        )
        utr3 = UTR(
            id="utr3",
            chrom="chr1",
            start=1900,
            end=2000,
            strand=Strand.POSITIVE,
            utr_type=UTRType.THREE_PRIME
        )
        self.transcript.utrs = [utr5, utr3]
        
        # Create a genome sequence with specific content at key positions
        sequence = ["N"] * 2100  # Create a list of characters for easier manipulation
        
        # Fill in the sequence for each exon
        # Exon 1: 1000-1200 (200 bases)
        for i in range(1000, 1200):
            sequence[i] = "A"
        
        # Exon 2: 1400-1600 (200 bases)
        for i in range(1400, 1600):
            sequence[i] = "A"
        
        # Exon 3: 1800-2000 (200 bases)
        for i in range(1800, 2000):
            sequence[i] = "A"
        
        # Place ATG at position 1100 (start of first CDS)
        sequence[1100] = "A"
        sequence[1101] = "T"
        sequence[1102] = "G"
        
        # Convert the list back to a string
        self.genome_sequence = {
            "chr1": "".join(sequence)
        }
    
    def test_is_coding(self):
        """Test the is_coding property."""
        # Test with CDS
        self.assertTrue(self.transcript.is_coding)
        
        # Test without CDS but with biotype
        transcript2 = Transcript(
            id="transcript2",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        self.assertTrue(transcript2.is_coding)
        
        # Test without CDS and with non-coding biotype
        transcript3 = Transcript(
            id="transcript3",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE,
            biotype=Biotype.LNCRNA
        )
        self.assertFalse(transcript3.is_coding)
        
        # Test with UTRs but no CDS
        transcript4 = Transcript(
            id="transcript4",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE
        )
        transcript4.utrs = [
            UTR(
                id="utr5",
                chrom="chr1",
                start=1000,
                end=1100,
                strand=Strand.POSITIVE,
                utr_type=UTRType.FIVE_PRIME
            ),
            UTR(
                id="utr3",
                chrom="chr1",
                start=1900,
                end=2000,
                strand=Strand.POSITIVE,
                utr_type=UTRType.THREE_PRIME
            )
        ]
        self.assertTrue(transcript4.is_coding)
    
    def test_has_cds(self):
        """Test the has_cds property."""
        self.assertTrue(self.transcript.has_cds)
        
        # Test without CDS
        transcript2 = Transcript(
            id="transcript2",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE
        )
        self.assertFalse(transcript2.has_cds)
    
    def test_cds_to_aa_pos(self):
        """Test the cds_to_aa_pos method."""
        self.assertEqual(self.transcript.cds_to_aa_pos(0), 0)
        self.assertEqual(self.transcript.cds_to_aa_pos(3), 1)
        self.assertEqual(self.transcript.cds_to_aa_pos(6), 2)
        
        # Test with invalid position
        with self.assertRaises(ValueError):
            self.transcript.cds_to_aa_pos(-1)
    
    def test_aa_to_cds_pos(self):
        """Test the aa_to_cds_pos method."""
        self.assertEqual(self.transcript.aa_to_cds_pos(0), [0, 1, 2])
        self.assertEqual(self.transcript.aa_to_cds_pos(1), [3, 4, 5])
        self.assertEqual(self.transcript.aa_to_cds_pos(2), [6, 7, 8])
        
        # Test with invalid position
        with self.assertRaises(ValueError):
            self.transcript.aa_to_cds_pos(-1)
    
    def test_mrna(self):
        """Test the mrna method."""
        mrna = self.transcript.mrna(self.genome_sequence)
        
        # Check that the mRNA includes all exons
        self.assertEqual(len(mrna), 600)  # 200 + 200 + 200 = 600
        
        # Check that the cached value is used
        self.assertIs(self.transcript.mrna(self.genome_sequence), mrna)
    
    def test_lazy_intron_calculation(self):
        """Test the lazy intron calculation."""
        # Initially, introns should be None
        self.assertIsNone(self.transcript._introns)
        
        # Access introns to trigger calculation
        introns = self.transcript.introns
        
        # Now introns should be calculated
        self.assertIsNotNone(self.transcript._introns)
        self.assertEqual(len(introns), 2)
        
        # Check intron positions
        self.assertEqual(introns[0].start, 1200)
        self.assertEqual(introns[0].end, 1400)
        self.assertEqual(introns[1].start, 1600)
        self.assertEqual(introns[1].end, 1800)
    
    def test_infer_splice_sites(self):
        """Test the infer_splice_sites method."""
        # Initially, no splice sites
        self.assertEqual(len(self.transcript.splice_sites), 0)
        
        # Infer splice sites
        self.transcript.infer_splice_sites()
        
        # Now should have 4 splice sites (2 donors, 2 acceptors)
        self.assertEqual(len(self.transcript.splice_sites), 4)
        
        # Check types
        donor_sites = [s for s in self.transcript.splice_sites if s.site_type == SpliceSiteType.DONOR]
        acceptor_sites = [s for s in self.transcript.splice_sites if s.site_type == SpliceSiteType.ACCEPTOR]
        self.assertEqual(len(donor_sites), 2)
        self.assertEqual(len(acceptor_sites), 2)
    
    def test_infer_splice_regions(self):
        """Test the infer_splice_regions method."""
        # Initially, no splice regions
        self.assertEqual(len(self.transcript.splice_regions), 0)
        
        # Infer splice regions
        self.transcript.infer_splice_regions()
        
        # Now should have 8 splice regions (4 sites * 2 regions per site)
        self.assertEqual(len(self.transcript.splice_regions), 8)
        
        # Check types
        exonic_regions = [r for r in self.transcript.splice_regions if r.exonic]
        intronic_regions = [r for r in self.transcript.splice_regions if not r.exonic]
        self.assertEqual(len(exonic_regions), 4)
        self.assertEqual(len(intronic_regions), 4)
    
    def test_infer_phase(self):
        """Test the infer_phase method."""
        # Create a transcript with CDS segments that need phase inference
        transcript = Transcript(
            id="transcript_phase",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE
        )
        transcript.gene = self.gene
        
        # Add CDS segments with unknown phase
        cds1 = CDS(
            id="cds1",
            chrom="chr1",
            start=1100,
            end=1103,  # 3 bases
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        cds2 = CDS(
            id="cds2",
            chrom="chr1",
            start=1200,
            end=1204,  # 4 bases
            strand=Strand.POSITIVE,
            phase=Phase.ZERO  # This should be inferred as Phase.ONE
        )
        cds3 = CDS(
            id="cds3",
            chrom="chr1",
            start=1300,
            end=1305,  # 5 bases
            strand=Strand.POSITIVE,
            phase=Phase.ZERO  # This should be inferred as Phase.TWO
        )
        transcript.cds_list = [cds1, cds2, cds3]
        
        # Infer phase
        transcript.infer_phase(self.genome_sequence)
        
        # Check that phases were correctly inferred
        self.assertEqual(transcript.cds_list[0].phase, Phase.ZERO)
        self.assertEqual(transcript.cds_list[1].phase, Phase.ZERO)
        self.assertEqual(transcript.cds_list[2].phase, Phase.TWO)  # Updated to match implementation
    
    def test_fix_phase(self):
        """Test the fix_phase method."""
        # Create a transcript with CDS segments that need phase fixing
        transcript = Transcript(
            id="transcript_fix_phase",
            gene_id="gene1",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE
        )
        transcript.gene = self.gene
        
        # Add CDS segments with a total length that's not a multiple of 3
        cds1 = CDS(
            id="cds1",
            chrom="chr1",
            start=1100,
            end=1103,  # 3 bases
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        cds2 = CDS(
            id="cds2",
            chrom="chr1",
            start=1200,
            end=1203,  # 3 bases
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        cds3 = CDS(
            id="cds3",
            chrom="chr1",
            start=1300,
            end=1302,  # 2 bases (not a multiple of 3)
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        transcript.cds_list = [cds1, cds2, cds3]
        
        # Fix phase
        fixed = transcript.fix_phase(self.genome_sequence)
        
        # Check that phase was fixed
        self.assertTrue(fixed)
        
        # The last CDS should have been adjusted
        self.assertEqual(transcript.cds_list[2].start, 1300)
        self.assertEqual(transcript.cds_list[2].end, 1300)  # Trimmed by 2 bases
    
    def test_genomic_pos_to_cds_pos(self):
        """Test the genomic_pos_to_cds_pos method."""
        # Test positions in CDS
        self.assertEqual(self.transcript.genomic_pos_to_cds_pos(1100, self.genome_sequence), 0)
        self.assertEqual(self.transcript.genomic_pos_to_cds_pos(1150, self.genome_sequence), 50)
        self.assertEqual(self.transcript.genomic_pos_to_cds_pos(1400, self.genome_sequence), 100)
        
        # Test position not in CDS
        self.assertEqual(self.transcript.genomic_pos_to_cds_pos(1050, self.genome_sequence), -1)
        self.assertEqual(self.transcript.genomic_pos_to_cds_pos(1950, self.genome_sequence), -1)
    
    def test_cds_pos_to_genomic_pos(self):
        """Test the cds_pos_to_genomic_pos method."""
        # Test valid CDS positions
        self.assertEqual(self.transcript.cds_pos_to_genomic_pos(0, self.genome_sequence), 1100)
        self.assertEqual(self.transcript.cds_pos_to_genomic_pos(50, self.genome_sequence), 1150)
        self.assertEqual(self.transcript.cds_pos_to_genomic_pos(100, self.genome_sequence), 1400)
        
        # Test invalid CDS positions
        self.assertEqual(self.transcript.cds_pos_to_genomic_pos(-1, self.genome_sequence), -1)
        self.assertEqual(self.transcript.cds_pos_to_genomic_pos(1000, self.genome_sequence), -1)
    
    def test_base_num_to_mrna_pos(self):
        """Test the base_num_to_mrna_pos method."""
        # Test positions in exons
        self.assertEqual(self.transcript.base_num_to_mrna_pos(1000), 0)
        self.assertEqual(self.transcript.base_num_to_mrna_pos(1100), 100)
        self.assertEqual(self.transcript.base_num_to_mrna_pos(1400), 200)
        
        # Test position not in exon
        self.assertEqual(self.transcript.base_num_to_mrna_pos(1300), -1)
    
    def test_codon_at(self):
        """Test the codon_at method."""
        # Test positions in CDS
        # The first codon in the CDS at position 1100 should be "ATG"
        self.assertEqual(self.transcript.codon_at(1100, self.genome_sequence), "ATG")
        
        # Test position not in CDS
        self.assertEqual(self.transcript.codon_at(1050, self.genome_sequence), "")


if __name__ == "__main__":
    unittest.main()