"""Tests for genomic annotation models."""

import unittest

from pygnome.genomics import (
    Biotype, CDS, Chromosome, Exon, Gene, Genome, Intron, 
    Phase, SpliceSite, SpliceSiteType, Strand, Transcript, UTR, UTRType
)


class TestGenomicModels(unittest.TestCase):
    """Test cases for genomic annotation models."""
    
    def test_basic_gene_model(self):
        """Test creating a basic gene model."""
        # Create a gene
        gene = Gene(
            id="ENSG00000139618",
            name="BRCA2",
            chrom="13",
            start=32315086,
            end=32400266,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        
        # Create a transcript
        transcript = Transcript(
            id="ENST00000380152",
            gene_id=gene.id,
            chrom=gene.chrom,
            start=32315086,
            end=32400266,
            strand=gene.strand,
            biotype=gene.biotype
        )
        
        # Add exons
        exon1 = Exon(
            id="ENSE00001484009",
            chrom=gene.chrom,
            start=32315086,
            end=32315358,
            strand=gene.strand
        )
        
        exon2 = Exon(
            id="ENSE00003659301",
            chrom=gene.chrom,
            start=32316422,
            end=32316527,
            strand=gene.strand
        )
        
        # Add CDS
        cds1 = CDS(
            id="CDS1",
            chrom=gene.chrom,
            start=32316461,
            end=32316527,
            strand=gene.strand,
            phase=Phase.ZERO
        )
        
        # Add UTR
        utr = UTR(
            id="UTR1",
            chrom=gene.chrom,
            start=32315086,
            end=32316460,
            strand=gene.strand,
            utr_type=UTRType.FIVE_PRIME
        )
        
        # Add to transcript
        transcript.exons = [exon1, exon2]
        transcript.cds_list = [cds1]
        transcript.utrs = [utr]
        
        # Add transcript to gene
        gene.transcripts = [transcript]
        
        # Verify properties
        self.assertEqual(gene.name, "BRCA2")
        self.assertEqual(gene.biotype, Biotype.PROTEIN_CODING)
        self.assertTrue(gene.is_coding)
        self.assertEqual(len(gene.transcripts), 1)
        self.assertEqual(len(transcript.exons), 2)
        self.assertEqual(transcript.coding_length, 66)  # 32316527 - 32316461
        self.assertEqual(len(transcript.five_prime_utrs), 1)
        self.assertEqual(len(transcript.three_prime_utrs), 0)
    
    def test_chromosome_and_genome(self):
        """Test creating chromosomes and genome."""
        # Create a simple gene
        gene = Gene(
            id="ENSG00000139618",
            name="BRCA2",
            chrom="13",
            start=32315086,
            end=32400266,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        
        # Create a chromosome
        chrom = Chromosome(name="13", length=114364328)
        chrom.add_gene(gene)
        
        # Create a genome
        genome = Genome(name="GRCh38", species="Homo sapiens")
        genome.add_chromosome(chrom)
        
        # Verify
        self.assertEqual(len(chrom.genes), 1)
        self.assertEqual(len(genome.chromosomes), 1)
        self.assertEqual(len(genome.genes), 1)
        self.assertEqual(genome.get_gene("ENSG00000139618"), gene)
        self.assertEqual(genome.get_chromosome("13"), chrom)
        
        # Test iteration
        chroms = list(genome)
        self.assertEqual(len(chroms), 1)
        self.assertEqual(chroms[0], chrom)
        
        genes = list(chrom)
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0], gene)


if __name__ == "__main__":
    unittest.main()