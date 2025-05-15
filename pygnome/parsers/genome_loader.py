"""
Genome Loader for building complete genomes from annotation and sequence files.
"""
from pathlib import Path
from tqdm import tqdm

from ..genomics.genome import Genome
from ..genomics.chromosome import Chromosome
from ..genomics.gene import Gene
from .gff.genome_loader import GffGenomeLoader
from .fasta.fasta_parser import FastaParser
from ..sequences.dna_string import DnaString


class GenomeLoader:
    """
    Class for loading complete genomes from annotation and sequence files.
    
    This class combines sequence data from FASTA files with annotation data
    from GFF/GTF files to build a complete Genome object with chromosomes,
    genes, transcripts, exons, and other genomic features.
    """
    
    def __init__(self, genome_name: str = "genome", species: str = None, verbose: bool = False):
        """
        Initialize the genome loader.
        
        Args:
            genome_name: Name of the genome
            species: Species name
            verbose: Whether to print progress information during loading
        """
        self.genome_name = genome_name
        self.species = species
        self.verbose = verbose
        self.genome = Genome(name=genome_name, species=species)
        self.feature_loader = GffGenomeLoader()
    
    def load(self,
             annotation_file: Path,
             sequence_file: Path) -> Genome:
        """
        Load a genome from annotation and sequence files.
        
        Args:
            annotation_file: Path to the GFF/GTF annotation file
            sequence_file: Path to the FASTA sequence file
            
        Returns:
            Loaded Genome object
        """
        # Load chromosome sequences
        if self.verbose:
            print(f"Loading chromosome sequences from {sequence_file}")
        chromosomes = self.load_sequences(sequence_file)
        
        # Load genomic features
        if self.verbose:
            print(f"Loading genomic features from {annotation_file}")
        self.feature_loader.load(annotation_file, self.verbose)
        
        # Add features to chromosomes and add chromosomes to genome
        for chrom in chromosomes.values():
            chrom_genes = [gene for gene in self.genes_by_id.values() if gene.chrom == chrom.name]
            chrom.add_genes(chrom_genes)
            self.genome.add_chromosome(chrom)
            
        if self.verbose:
            print(f"Genome loading complete: {len(self.genome.chromosomes)} chromosomes, {len(self.genome.genes)} genes")
        
        return self.genome

    def load_sequences(self, sequence_file: Path) -> dict[str, Chromosome]:
        """
        Load chromosome sequences from a FASTA file.
        
        Args:
            sequence_file: Path to the FASTA file
            
        Returns:
            Dictionary of chromosome names to Chromosome objects
        """
        # Load sequences using the FASTA parser
        dna_strings = FastaParser.parse_as_dna_strings(sequence_file)
        if not dna_strings:
            raise ValueError(f"No sequences found in {sequence_file}")
        
        # Create chromosomes from sequences
        chromosomes = {}
        
        if self.verbose and len(dna_strings) > 1:
            # Use progress bar for multiple chromosomes
            chrom_iter = tqdm(dna_strings.items(), desc="Loading chromosomes")
        else:
            chrom_iter = dna_strings.items()
            
        for chrom_name, dna_sequence in chrom_iter:
            chromosome = Chromosome(name=chrom_name, sequence=dna_sequence)
            chromosomes[chrom_name] = chromosome
        
        return chromosomes
        
    def __str__(self) -> str:
        """
        Return a string representation of the genome loader.
        
        This method provides a summary of the loaded genome without
        attempting to print the entire genome data, which could be
        several gigabytes in size.
        """
        genome_str = f"GenomeLoader({self.genome_name}"
        if self.species:
            genome_str += f", {self.species}"
        genome_str += ")"
        
        if not self.genome.chromosomes:
            return f"{genome_str} - No chromosomes loaded"
        
        chrom_count = len(self.genome.chromosomes)
        gene_count = len(self.genome.genes)
        
        # Calculate total sequence length
        total_length = 0
        for chrom in self.genome.chromosomes.values():
            if chrom.length is not None:
                total_length += chrom.length
        
        # Format total length in a human-readable way
        if total_length > 1_000_000_000:
            length_str = f"{total_length / 1_000_000_000:.2f} Gb"
        elif total_length > 1_000_000:
            length_str = f"{total_length / 1_000_000:.2f} Mb"
        elif total_length > 1_000:
            length_str = f"{total_length / 1_000:.2f} kb"
        else:
            length_str = f"{total_length} bp"
        
        return f"{genome_str} - {chrom_count} chromosomes, {gene_count} genes, {length_str}"