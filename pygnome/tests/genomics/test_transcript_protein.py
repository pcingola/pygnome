"""Test the Transcript protein sequence functionality."""

import unittest

from pygnome.genomics.biotype import Biotype
from pygnome.genomics.cds import CDS
from pygnome.genomics.chromosome import Chromosome
from pygnome.genomics.codon_table import CodonTable, CodonTableType
from pygnome.genomics.gene import Gene
from pygnome.genomics.genome import Genome
from pygnome.genomics.phase import Phase
from pygnome.genomics.strand import Strand
from pygnome.genomics.transcript import Transcript


class TestTranscriptProtein(unittest.TestCase):
    """Test the Transcript protein sequence functionality."""
    
    def test_protein_with_genome_codon_table(self):
        """Test that the transcript uses the genome's codon table."""
        # Create a genome sequence with ATA codon (codes for I in standard, M in vertebrate mitochondrial)
        # Make sure the sequence is long enough to cover the CDS positions (100-109)
        genome_sequence = {"chr1": "." * 100 + "ATAATAATA"}  # Three ATA codons at positions 100-109
        
        # Test with standard codon table
        standard_genome = Genome(
            name="Standard Genome",
            species="Test Species",
            codon_table_type=CodonTableType.STANDARD
        )
        
        standard_chrom = Chromosome(name="chr1")
        standard_genome.add_chromosome(standard_chrom)
        
        standard_gene = Gene(
            id="gene1",
            chrom="chr1",
            start=100,
            end=200,
            strand=Strand.POSITIVE,
            name="Gene1"
        )
        standard_chrom.add_gene(standard_gene)
        
        standard_transcript = Transcript(
            id="transcript1",
            gene_id="gene1",
            chrom="chr1",
            start=100,
            end=109,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        
        standard_cds = CDS(
            id="cds1",
            chrom="chr1",
            start=100,
            end=109,  # 9 bases = 3 codons
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        standard_transcript.cds_list = [standard_cds]
        standard_gene.transcripts = [standard_transcript]
        standard_transcript.gene = standard_gene
        
        # In standard code, ATA codes for I, so should be "III"
        protein_standard = standard_transcript.protein(genome_sequence)
        self.assertEqual(protein_standard, "III")
        
        # Test with vertebrate mitochondrial codon table
        mito_genome = Genome(
            name="Mitochondrial Genome",
            species="Test Species",
            codon_table_type=CodonTableType.VERTEBRATE_MITOCHONDRIAL
        )
        
        mito_chrom = Chromosome(name="chr1")
        mito_genome.add_chromosome(mito_chrom)
        
        mito_gene = Gene(
            id="gene1",
            chrom="chr1",
            start=100,
            end=200,
            strand=Strand.POSITIVE,
            name="Gene1"
        )
        mito_chrom.add_gene(mito_gene)
        
        mito_transcript = Transcript(
            id="transcript1",
            gene_id="gene1",
            chrom="chr1",
            start=100,
            end=109,
            strand=Strand.POSITIVE,
            biotype=Biotype.PROTEIN_CODING
        )
        
        mito_cds = CDS(
            id="cds1",
            chrom="chr1",
            start=100,
            end=109,  # 9 bases = 3 codons
            strand=Strand.POSITIVE,
            phase=Phase.ZERO
        )
        mito_transcript.cds_list = [mito_cds]
        mito_gene.transcripts = [mito_transcript]
        mito_transcript.gene = mito_gene
        
        # In vertebrate mitochondrial code, ATA codes for M, so should be "MMM"
        protein_mito = mito_transcript.protein(genome_sequence)
        self.assertEqual(protein_mito, "MMM")