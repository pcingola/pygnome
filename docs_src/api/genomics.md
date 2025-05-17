# Genomics

The `genomics` module provides the core data models for representing genomic features. It includes classes for representing genomes, chromosomes, genes, transcripts, and other genomic features.

## Overview

The genomics module is built around a hierarchy of classes that represent different levels of genomic organization:

1. `Genome`: Top-level container for chromosomes and genes
2. `Chromosome`: Contains genes and sequence data
3. `Gene`: Represents a gene with transcripts
4. `Transcript`: Represents a transcript with exons, introns, UTRs, and CDS
5. `Exon`, `Intron`, `UTR`, `CDS`: Represent specific genomic features

All genomic features inherit from the `GenomicFeature` base class, which provides common functionality for working with genomic coordinates.

## Core Classes

### GenomicFeature

```python
@dataclass
class GenomicFeature:
```

Base class for all genomic features. Uses 0-based, half-open interval coordinates [start, end) where start is included but end is excluded. Length is calculated as (end - start).

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | `str` | Unique identifier for the feature |
| `chrom` | `str` | Chromosome name |
| `start` | `int` | Start position (0-based, inclusive) |
| `end` | `int` | End position (0-based, exclusive) |
| `strand` | `Strand` | Strand orientation (positive/negative) |

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `length` | `int` | Length of the feature (end - start) |

#### Methods

| Method | Description |
|--------|-------------|
| `intersects_point(position: int) -> bool` | Check if the feature intersects with a specific point |
| `intersects_interval(start: int, end: int) -> bool` | Check if the feature intersects with a specific interval |
| `intersects(other: GenomicFeature) -> bool` | Check if the feature intersects with another feature |
| `contains(other: GenomicFeature) -> bool` | Check if the feature contains another feature |
| `distance(position: int) -> int` | Calculate the distance from a point to the feature |

### Genome

```python
class Genome:
```

A genome containing chromosomes and genes.

#### Constructor

```python
def __init__(self, name: str, species: str | None = None, codon_table_type: CodonTableType = CodonTableType.STANDARD):
```

- `name`: Name of the genome
- `species`: Species name (optional)
- `codon_table_type`: Type of codon table to use for translation (default: `CodonTableType.STANDARD`)

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Name of the genome |
| `species` | `str \| None` | Species name |
| `codon_table_type` | `CodonTableType` | Type of codon table to use for translation |
| `chromosomes` | `dict[str, Chromosome]` | Dictionary of chromosomes by name |
| `genes` | `dict[str, Gene]` | Dictionary of genes by ID |

#### Methods

| Method | Description |
|--------|-------------|
| `add_chromosome(chromosome: Chromosome) -> None` | Add a chromosome to the genome |
| `get_chromosome(chrom_name: str) -> Chromosome \| None` | Get a chromosome by name |
| `get_gene(gene_id: str) -> Gene \| None` | Get a gene by ID |
| `get_gene_by_name(gene_name: str) -> list[Gene]` | Get genes by name |
| `get_genes_by_biotype(biotype: str) -> list[Gene]` | Get genes by biotype |

### Chromosome

```python
class Chromosome:
```

A chromosome containing genes and optionally a DNA sequence.

#### Constructor

```python
def __init__(self, name: str, length: int | None = None, sequence: DnaString | None = None):
```

- `name`: Name of the chromosome
- `length`: Length of the chromosome in base pairs (optional)
- `sequence`: DNA sequence of the chromosome (optional)

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Name of the chromosome |
| `sequence` | `DnaString \| None` | DNA sequence of the chromosome |
| `genes` | `dict[str, Gene]` | Dictionary of genes by ID |
| `genome` | `Genome \| None` | Reference to the parent genome |

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `length` | `int \| None` | Length of the chromosome in base pairs |

#### Methods

| Method | Description |
|--------|-------------|
| `add_gene(gene: Gene) -> None` | Add a gene to the chromosome |
| `add_genes(genes: list[Gene]) -> None` | Add multiple genes to the chromosome |
| `get_gene(gene_id: str) -> Gene \| None` | Get a gene by ID |
| `get_genes_in_region(start: int, end: int) -> list[Gene]` | Get all genes that overlap a genomic region |

### Gene

```python
@dataclass
class Gene(GenomicFeature):
```

A gene, which may have multiple transcripts.

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str \| None` | Gene name |
| `biotype` | `Biotype \| None` | Gene biotype |
| `transcripts` | `list[Transcript]` | List of transcripts |
| `chromosome` | `Chromosome \| None` | Reference to the parent chromosome |

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `canonical_transcript` | `Transcript \| None` | Return the canonical transcript, if defined |
| `is_coding` | `bool` | Return True if any transcript has coding sequence |

#### Methods

| Method | Description |
|--------|-------------|
| `add_transcript(transcript: Transcript) -> None` | Add a transcript to the gene |

### Transcript

```python
@dataclass
class Transcript(GenomicFeature):
```

A transcript of a gene, composed of exons, introns, UTRs, and CDS segments.

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `gene_id` | `str` | ID of the parent gene |
| `biotype` | `Biotype \| None` | Transcript biotype |
| `exons` | `list[Exon]` | List of exons |
| `cds_list` | `list[CDS]` | List of CDS segments |
| `utrs` | `list[UTR]` | List of UTRs |
| `introns` | `list[Intron]` | List of introns |
| `splice_sites` | `list[SpliceSite]` | List of splice sites |
| `gene` | `Gene \| None` | Reference to the parent gene |

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `coding_length` | `int` | Total length of coding sequence |
| `exonic_length` | `int` | Total length of exons |
| `is_coding` | `bool` | Return True if the transcript has coding sequence |
| `five_prime_utrs` | `list[UTR]` | Return all 5' UTRs |
| `three_prime_utrs` | `list[UTR]` | Return all 3' UTRs |

#### Methods

| Method | Description |
|--------|-------------|
| `cds(genome_sequence: dict[str, str]) -> str` | Get the coding DNA sequence for this transcript |
| `protein(genome_sequence: dict[str, str]) -> str` | Get the protein sequence for this transcript |

## Enumerations

### Strand

```python
class Strand(str, Enum):
```

Enumeration for strand orientation.

| Value | Description |
|-------|-------------|
| `POSITIVE` | Positive strand (+) |
| `NEGATIVE` | Negative strand (-) |
| `UNKNOWN` | Unknown strand (.) |

### Biotype

```python
class Biotype(str, Enum):
```

Enumeration for gene/transcript biotypes.

| Value | Description |
|-------|-------------|
| `PROTEIN_CODING` | Protein-coding gene/transcript |
| `NONCODING` | Non-coding gene/transcript |
| `PSEUDOGENE` | Pseudogene |
| `UNKNOWN` | Unknown biotype |

### Phase

```python
class Phase(int, Enum):
```

Enumeration for CDS phase.

| Value | Description |
|-------|-------------|
| `ZERO` | Phase 0 (0) |
| `ONE` | Phase 1 (1) |
| `TWO` | Phase 2 (2) |
| `UNKNOWN` | Unknown phase (-1) |

## Usage Examples

### Creating a Genome

```python
from pygnome.genomics.genome import Genome
from pygnome.genomics.chromosome import Chromosome
from pygnome.genomics.gene import Gene
from pygnome.genomics.transcript import Transcript
from pygnome.genomics.exon import Exon
from pygnome.genomics.strand import Strand

# Create a genome
genome = Genome(name="MyGenome", species="Example species")

# Create a chromosome
chromosome = Chromosome(name="chr1", length=248956422)

# Add the chromosome to the genome
genome.add_chromosome(chromosome)

# Create a gene
gene = Gene(
    id="GENE001",
    name="Example Gene",
    chrom="chr1",
    start=1000,
    end=5000,
    strand=Strand.POSITIVE
)

# Create a transcript
transcript = Transcript(
    id="TRANSCRIPT001",
    gene_id="GENE001",
    chrom="chr1",
    start=1000,
    end=5000,
    strand=Strand.POSITIVE
)

# Create exons
exon1 = Exon(
    id="EXON001",
    chrom="chr1",
    start=1000,
    end=2000,
    strand=Strand.POSITIVE
)

exon2 = Exon(
    id="EXON002",
    chrom="chr1",
    start=3000,
    end=5000,
    strand=Strand.POSITIVE
)

# Add exons to transcript
transcript.exons = [exon1, exon2]

# Add transcript to gene
gene.add_transcript(transcript)

# Add gene to chromosome
chromosome.add_gene(gene)

# Now the gene is also accessible from the genome
print(genome.genes["GENE001"])
```

### Working with Genomic Features

```python
from pygnome.genomics.gene import Gene
from pygnome.genomics.strand import Strand

# Create two genes
gene1 = Gene(
    id="GENE001",
    name="Gene 1",
    chrom="chr1",
    start=1000,
    end=5000,
    strand=Strand.POSITIVE
)

gene2 = Gene(
    id="GENE002",
    name="Gene 2",
    chrom="chr1",
    start=4000,
    end=7000,
    strand=Strand.NEGATIVE
)

# Check if the genes intersect
if gene1.intersects(gene2):
    print(f"{gene1.id} intersects with {gene2.id}")

# Check if a position is within a gene
position = 3000
if gene1.intersects_point(position):
    print(f"Position {position} is within {gene1.id}")

# Calculate the distance from a position to a gene
position = 6000
distance = gene1.distance(position)
print(f"Distance from position {position} to {gene1.id}: {distance}")
```