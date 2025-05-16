# PyGnome: Python Library for Genome Annotations

PyGnome is a Python library for working with genomic annotations, similar to tools like SnpEff/SnpSift but as a Python library.

## Project Structure

```
pygnome/
├── genomics/               # Genomic feature models
│   ├── gene.py             # Gene class
│   ├── transcript.py       # Transcript class
│   ├── exon.py             # Exon class
│   └── feature_store/      # Feature storage implementations
├── parsers/                # File format parsers
│   ├── fasta/              # FASTA/FASTQ parsers
│   ├── gff/                # GFF/GTF parsers
│   ├── msi/                # MSI-sites parsers
│   └── vcf/                # VCF parsers
├── sequences/              # Sequence handling
│   ├── base_sequence.py    # Base sequence class
│   ├── dna_string.py       # DNA sequence implementation
│   └── rna_string.py       # RNA sequence implementation
└── tests/                  # Test suite
```

## Components

### Genomic Features
The library models genomic features with a hierarchical structure:

- `GenomicFeature`: Base class for all genomic features with coordinates (chrom, start, end, strand)
- `Gene`: Contains multiple transcripts and gene-level properties
- `Transcript`: Represents an RNA transcript with exons, introns, UTRs, and coding sequences
- `Exon`, `Intron`, `CDS`, `UTR`, `SpliceSite`: Specific genomic features with specialized properties
- `Variant`: Base class for genomic variants with reference and alternate alleles
  - `SNP`: Single Nucleotide Polymorphism (single base change)
  - `Insertion`: Insertion of one or more bases
  - `Deletion`: Deletion of one or more bases
  - `Duplication`: Duplication of a genomic segment
  - `Inversion`: Inversion of a genomic segment
  - `Translocation`: Movement of a genomic segment to a different location
  - `ComplexVariant`: Complex variants that don't fit other categories
- `Strand`, `Phase`, `Biotype`, `UTRType`: Enumerations for genomic feature properties

### Feature Storage
The library provides efficient data structures for querying genomic features:

- `FeatureStore`: Interface for storing and querying genomic features
- `GenomicFeatureStore`: Main implementation with multiple storage backends:
  - `IntervalTreeStore`: Uses interval trees for efficient range queries
  - `BinnedGenomicStore`: Divides genome into bins for faster lookups
  - `BruteForceFeatureStore`: Simple linear search implementation
  - `MsiChromosomeStore`: Specialized store for MSI (Microsatellite Instability) sites
- Feature stores support serialization via pickle for quick saving and loading of genomic data

### Parsers
Parsers for common genomic file formats:

- GFF/GTF Parsers: Parse Gene Feature Format files (GFF2, GFF3, GTF)
- FASTA/FASTQ Parsers: Parse sequence files
- MSI-sites Parser: Parse Microsatellite Instability sites files
- VCF Parsers: Parse Variant Call Format files
  - `VcfReader`: Read and parse VCF files
  - `VcfHeader`: Parse and represent VCF headers
  - `VcfRecord`: Represent individual VCF records
  - `VariantFactory`: Create Variant objects from VCF records

### Sequence Handling
Efficient representation of nucleotide sequences:

- `BaseSequence`: Abstract base class for 2-bit encoded nucleotide sequences
- `DnaString`: DNA-specific sequence implementation (A,C,G,T)
- `RnaString`: RNA-specific sequence implementation (A,C,G,U)

The library uses Pydantic for data validation and provides a comprehensive test suite for all components. It's designed to efficiently handle genomic data with a focus on annotation and feature querying.

