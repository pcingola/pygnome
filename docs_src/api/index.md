# API Reference

This section provides detailed documentation for the PyGnome API. PyGnome is organized into several core modules, each focusing on a specific aspect of genomic data handling.

## Core Modules

### [Genomics](genomics.md)

The `genomics` module provides the core data models for representing genomic features:

- `GenomicFeature`: Base class for all genomic features
- `Genome`: Top-level container for chromosomes and genes
- `Chromosome`: Contains genes and sequence data
- `Gene`: Represents a gene with transcripts
- `Transcript`: Represents a transcript with exons, introns, UTRs, and CDS
- `Exon`, `Intron`, `UTR`, `CDS`: Represent specific genomic features
- `Strand`: Enumeration for strand orientation (positive/negative)
- `Biotype`: Enumeration for gene/transcript biotypes
- `Phase`: Enumeration for CDS phase

### [Feature Store](feature-store.md)

The `feature_store` module provides efficient storage and retrieval of genomic features:

- `GenomicFeatureStore`: Main interface for storing and querying features
- `IntervalTreeStore`: Uses interval trees for efficient range queries
- `BinnedGenomicStore`: Uses binning for memory-efficient storage
- `BruteForceFeatureStore`: Simple implementation for testing
- `MsiChromosomeStore`: Specialized for microsatellite instability sites

### [Sequences](sequences.md)

The `sequences` module provides memory-efficient representations of DNA and RNA sequences:

- `BaseSequence`: Abstract base class for nucleotide sequences
- `DnaString`: Memory-efficient 2-bit representation of DNA sequences
- `RnaString`: RNA sequence representation
- `DnaStringArray`: Efficient storage for multiple DNA sequences

### [Parsers](parsers.md)

The `parsers` module includes parsers for common genomic file formats:

- `GenomeLoader`: Loads genomes from annotation and sequence files
- `FastaParser`: Parses FASTA files
- `FastqParser`: Parses FASTQ files
- `GffParser`, `Gff3Parser`, `GtfParser`: Parse GFF/GTF annotation files
- `VcfReader`: Parses VCF variant files
- `MsiSitesReader`: Parses microsatellite instability sites

## Module Relationships

The modules in PyGnome are designed to work together:

1. **Parsers** read genomic data from files
2. **Genomics** models represent the parsed data as objects
3. **Feature Stores** provide efficient storage and retrieval of genomic features
4. **Sequences** provide memory-efficient representation of DNA/RNA sequences

## Common Usage Patterns

### Loading and Querying a Genome

```python
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore

# Load a genome
loader = GenomeLoader(genome_name="GRCh38", species="Homo sapiens")
genome = loader.load(
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz")
)

# Create a feature store for efficient querying
store = GenomicFeatureStore()
with store:
    for gene in genome.genes.values():
        store.add(gene)
        for transcript in gene.transcripts:
            store.add(transcript)
            for exon in transcript.exons:
                store.add(exon)

# Query features
features = store.get_by_interval("chr1", 1000000, 2000000)
```

### Working with Sequences

```python
from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString

# Create a DNA sequence
dna = DnaString("ATGCATGCATGC")

# Get a subsequence
subseq = dna[3:9]

# Complement and reverse complement
comp = dna.complement()
rev_comp = dna.reverse_complement()

# Transcribe DNA to RNA
rna = dna.transcribe()

# Translate RNA to protein
protein = rna.translate()
```

## Error Handling

PyGnome uses Python's built-in exception handling. Common exceptions include:

- `ValueError`: Raised when a function receives an argument of the correct type but an inappropriate value
- `TypeError`: Raised when an operation or function is applied to an object of inappropriate type
- `IndexError`: Raised when a sequence subscript is out of range
- `FileNotFoundError`: Raised when a file or directory is requested but doesn't exist

Example of error handling:

```python
from pathlib import Path
from pygnome.parsers.fasta.fasta_parser import FastaParser

try:
    parser = FastaParser(Path("path/to/nonexistent_file.fa"))
    records = parser.load()
except FileNotFoundError as e:
    print(f"Error: {e}")
    # Handle the error appropriately