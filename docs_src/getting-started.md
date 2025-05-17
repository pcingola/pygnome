# Getting Started with PyGnome

PyGnome is a Python library for working with genomic annotations and sequences. It provides efficient data structures and parsers for common genomic file formats, making it easy to work with genomic data in Python.

## Installation

PyGnome can be installed using pip:

```bash
pip install pygnome
```

## Quick Start

Here's a simple example to get you started with PyGnome:

```python
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader

# Load a genome from GTF and FASTA files
loader = GenomeLoader(genome_name="GRCh38", species="Homo sapiens")
genome = loader.load(
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz")
)

# Access genomic features
for gene in genome.genes.values():
    print(f"Gene: {gene.id} ({gene.name}) - {gene.chrom}:{gene.start}-{gene.end}")
    
    for transcript in gene.transcripts:
        print(f"  Transcript: {transcript.id} - Exons: {len(transcript.exons)}")
```

## Core Components

PyGnome consists of several core components:

### Genomic Models

The genomic models provide a comprehensive object-oriented representation of genomic features:

- `Genome`: Top-level container for chromosomes and genes
- `Chromosome`: Contains genes and sequence data
- `Gene`: Represents a gene with transcripts
- `Transcript`: Represents a transcript with exons, introns, UTRs, and CDS
- `Exon`, `Intron`, `UTR`, `CDS`: Represent specific genomic features

### Feature Stores

Feature stores provide efficient storage and retrieval of genomic features:

- `GenomicFeatureStore`: Main interface for storing and querying features
- Multiple implementations with different performance characteristics:
  - `IntervalTreeStore`: Uses interval trees for efficient range queries
  - `BinnedGenomicStore`: Uses binning for memory-efficient storage
  - `BruteForceFeatureStore`: Simple implementation for testing
  - `MsiChromosomeStore`: Specialized for microsatellite instability sites

### Sequence Handling

PyGnome provides memory-efficient representations of DNA and RNA sequences:

- `DnaString`: Memory-efficient 2-bit representation of DNA sequences
- `RnaString`: RNA sequence representation
- `DnaStringArray`: Efficient storage for multiple DNA sequences

### Parsers

PyGnome includes parsers for common genomic file formats:

- `FastaParser`: Parses FASTA files
- `FastqParser`: Parses FASTQ files
- `GffParser`, `Gff3Parser`, `GtfParser`: Parse GFF/GTF annotation files
- `VcfReader`: Parses VCF variant files
- `MsiSitesReader`: Parses microsatellite instability sites

## Next Steps

- Check out the [User Guide](user-guide/basic-usage.md) for more detailed usage examples
- Explore the [API Reference](api/index.md) for detailed documentation of all classes and methods
- Learn about the [Format Specifications](format_specifications/fasta.md) for supported file formats