"""Genome class for genomic annotations."""

import re
from .chromosome import Chromosome
from .chromosome_category import ChromosomeCategory
from .codon_table import CodonTableType
from .gene import Gene


class Genome:
    """A genome containing chromosomes and genes."""
    
    def __init__(self, name: str, species: str | None = None,
                codon_table_type: CodonTableType = CodonTableType.STANDARD):
        """Initialize a genome."""
        self.name = name
        self.species = species
        self.codon_table_type = codon_table_type
        self.chromosomes: dict[str, Chromosome] = {}
        self.genes: dict[str, Gene] = {}
        
    def add_chromosome(self, chromosome: Chromosome) -> None:
        """Add a chromosome to the genome."""
        self.chromosomes[chromosome.name] = chromosome
        
        # Set the genome reference in the chromosome
        chromosome.genome = self
        
        # Add all genes from the chromosome to the genome's gene index
        for gene_id, gene in chromosome.genes.items():
            self.genes[gene_id] = gene
    
    def get_chromosome(self, chrom_name: str) -> Chromosome | None:
        """Get a chromosome by name."""
        return self.chromosomes.get(chrom_name)
    
    def get_gene(self, gene_id: str) -> Gene | None:
        """Get a gene by ID."""
        return self.genes.get(gene_id)
    
    def get_gene_by_name(self, gene_name: str) -> list[Gene]:
        """Get genes by name."""
        return [gene for gene in self.genes.values() if gene.name == gene_name]
    
    def get_genes_by_biotype(self, biotype: str) -> list[Gene]:
        """Get genes by biotype."""
        return [gene for gene in self.genes.values() if gene.biotype == biotype]
    
    def __str__(self) -> str:
        """Return a string representation of the genome."""
        species_str = f", {self.species}" if self.species else ""
        codon_table_str = f", {self.codon_table_type}"
        return f"Genome({self.name}{species_str}{codon_table_str}, {len(self.chromosomes)} chromosomes, {len(self.genes)} genes)"
    
    def __iter__(self):
        """Iterate over chromosomes in standard order (1-22, X, Y, M, others)."""
        sorted_chroms = sorted(self.chromosomes.values(), 
                              key=lambda c: self._chromosome_sort_key(c.name))
        return iter(sorted_chroms)
    
    def _chromosome_sort_key(self, chrom_name: str) -> tuple:
        """Create a sort key for chromosome names."""
        # Extract number if present (e.g., "chr1" -> 1, "chrX" -> None)
        match = re.match(r'(?:chr)?(\d+)$', chrom_name, re.IGNORECASE)
        if match:
            return (ChromosomeCategory.NUMERIC, int(match.group(1)))
        
        # Sex chromosomes (X, Y)
        if chrom_name.upper() in ('X', 'CHRX'):
            return (ChromosomeCategory.SEX, 0)  # X comes before Y
        if chrom_name.upper() in ('Y', 'CHRY'):
            return (ChromosomeCategory.SEX, 1)  # Y comes after X
        
        # Mitochondrial chromosome (M, MT)
        if chrom_name.upper() in ('M', 'MT', 'CHRM', 'CHRMT'):
            return (ChromosomeCategory.MITOCHONDRIAL, 0)
        
        # Other chromosomes sorted alphabetically
        return (ChromosomeCategory.OTHER, chrom_name)