
# PyGnome: Python Library for Genome Annotations

<img src="./docs_src/images/pygnome_alpha.png" width="200">


PyGnome is a Python library for working with genomic annotations and sequences. It provides efficient data structures and parsers for common genomic file formats, making it easy to work with genomic data in Python. 

At its core, PyGnome offers *"Genomic feature stores"*, which are specialized data structures for efficient storage, indexing, and querying of genomic features based on their coordinates, solving the fundamental bioinformatics challenge of quickly locating genomic elements within large genomes.

Full documentation is available at [https://pcingola.github.io/pygnome](https://pcingola.github.io/pygnome)

## Features

- **Genomic Feature Models**: Comprehensive object models for genes, transcripts, exons, variants, and more
- **Efficient Feature Storage**: Multiple implementations for fast genomic feature queries
- **File Format Parsers**: Support for FASTA/FASTQ, GFF/GTF, VCF, and MSI formats
- **Sequence Handling**: Memory-efficient representations of DNA and RNA sequences

## Installation

```bash
pip install pygnome
```

## Quick Start

```python
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader

# Load a genome from GTF and FASTA files
loader = GenomeLoader(genome_name="GRCh38", species="Homo sapiens")
genome = loader.load(
    gtf_file=Path("path/to/annotations.gtf"),
    fasta_file=Path("path/to/genome.fa.gz")
)

# Access genomic features
for gene in genome.genes.values():
    print(f"Gene: {gene.id} ({gene.name}) - {gene.chrom}:{gene.start}-{gene.end}")
    
    for transcript in gene.transcripts:
        print(f"  Transcript: {transcript.id} - Exons: {len(transcript.exons)}")
```

## Usage Examples

### Parsing FASTA Files

```python
from pathlib import Path
from pygnome.parsers.fasta.fasta_parser import FastaParser

# Parse a FASTA file
parser = FastaParser(Path("path/to/sequences.fa"))
records = parser.load()

# Access sequences
for record in records:
    print(f"Sequence: {record.identifier}")
    print(f"Length: {len(record.sequence)}")
    
    # Convert to string if needed
    seq_str = str(record.sequence)
    print(f"First 10 bases: {seq_str[:10]}")

# Load as dictionary for quick access by identifier
sequences = FastaParser(Path("path/to/sequences.fa")).load_as_dict()
my_seq = sequences["chr1"].sequence
```

### Parsing GFF/GTF Files

```python
from pathlib import Path
from pygnome.parsers.gff.gff3_parser import Gff3Parser
from pygnome.parsers.gff.gtf_parser import GtfParser

# Parse a GFF3 file
gff_parser = Gff3Parser(Path("path/to/annotations.gff3"))
for record in gff_parser:
    print(f"{record.type}: {record.chrom}:{record.start}-{record.end}")
    print(f"Attributes: {record.attributes}")

# Parse a GTF file
gtf_parser = GtfParser(Path("path/to/annotations.gtf"))
for record in gtf_parser:
    if record.type == "gene":
        gene_id = record.attributes.get("gene_id")
        gene_name = record.attributes.get("gene_name")
        print(f"Gene: {gene_id} ({gene_name}) - {record.chrom}:{record.start}-{record.end}")
```

### Parsing VCF Files

```python
from pathlib import Path
from pygnome.parsers.vcf.vcf_reader import VcfReader

# Open a VCF file
with VcfReader(Path("path/to/variants.vcf")) as reader:
    # Get sample names
    samples = reader.get_samples()
    print(f"Samples: {samples}")
    
    # Iterate through records
    for record in reader:
        print(f"Record: {record.get_chrom()}:{record.get_pos()} {record.get_ref()}>{','.join(record.get_alt())}")
        
        # Create variant objects from the record using VariantFactory
        for variant in record:  # Uses VariantFactory internally
            print(f"Variant: {variant}")
            
        # Access genotypes
        genotypes = record.get_genotypes()
        for i, genotype in enumerate(genotypes):
            print(f"  {samples[i]}: {genotype}")
        
    # Query a specific region
    for record in reader.fetch("chr1", 1000000, 2000000):
        for variant in record:
            print(f"Region variant: {variant}")
```

### Parsing VCF Annotations (ANN Field)

```python
from pathlib import Path
from pygnome.parsers.vcf.vcf_reader import VcfReader
from pygnome.parsers.vcf.ann import AnnParser

# Open a VCF file
with VcfReader(Path("path/to/variants.vcf")) as reader:
    # Iterate through records
    for record in reader:
        # Parse ANN field if present
        ann_parser = AnnParser(record)
        
        # Iterate through annotations
        for annotation in ann_parser:
            print(f"Variant annotation: {annotation.allele} - {annotation.annotation}")
            print(f"  Impact: {annotation.putative_impact}")
            
            if annotation.gene_name:
                print(f"  Gene: {annotation.gene_name}")
                
            if annotation.feature_type and annotation.feature_id:
                print(f"  Feature: {annotation.feature_type.value} {annotation.feature_id}")
                
            if annotation.hgvs_c:
                print(f"  HGVS.c: {annotation.hgvs_c}")
                
            if annotation.hgvs_p:
                print(f"  HGVS.p: {annotation.hgvs_p}")
```

### Using Genomic Feature Stores

Genomic feature stores are one of the core solutions in PyGnome, providing specialized data structures for efficient storage and retrieval of genomic features based on their genomic coordinates. They solve the fundamental bioinformatics challenge of quickly locating genomic elements within large genomes, allowing you to:

- Find all features at a specific position
- Find all features that overlap with a given range
- Find the nearest feature to a specific position

PyGnome offers multiple implementations with different performance characteristics to suit various use cases:

- **IntervalTreeStore** (default): Uses interval trees for efficient range queries
- **BinnedGenomicStore**: Uses binning for memory-efficient storage
- **BruteForceFeatureStore**: Simple implementation for testing
- **MsiChromosomeStore**: Specialized for microsatellite instability sites

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.genomics.gene import Gene
from pathlib import Path

# Create a feature store using interval trees (default)
store = GenomicFeatureStore()

# Or choose a different implementation
binned_store = GenomicFeatureStore(store_type=StoreType.BINNED, bin_size=100000)
brute_force_store = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)

# Add features to the store
with store:  # Use context manager to ensure proper indexing
    for gene in genome.genes.values():
        store.add(gene)
        
        # Add transcripts and other features
        for transcript in gene.transcripts:
            store.add(transcript)
            for exon in transcript.exons:
                store.add(exon)

# Query features
features_at_position = store.get_by_position("chr1", 1000000)
features_in_range = store.get_by_interval("chr1", 1000000, 2000000)
nearest_feature = store.get_nearest("chr1", 1500000)

# Save and load the store
store.save(Path("path/to/store.pkl"))
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
```

### Working with DNA/RNA Sequences

```python
from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString

# Create a DNA sequence
dna = DnaString("ATGCATGCATGC")
print(f"Length: {len(dna)}")
print(f"GC content: {dna.gc_content()}")

# Get a subsequence
subseq = dna[3:9]  # Returns a new DnaString

# Complement and reverse complement
comp = dna.complement()
rev_comp = dna.reverse_complement()

# Transcribe DNA to RNA
rna = dna.transcribe()  # Returns an RnaString

# Create an RNA sequence
rna = RnaString("AUGCAUGCAUGC")

# Translate RNA to protein
protein = rna.translate()
print(f"Protein: {protein}")
```

## Advanced Usage

### Loading a Complete Genome

```python
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader

# Create a genome loader
loader = GenomeLoader(
    genome_name="GRCh38",
    species="Homo sapiens",
    verbose=True  # Print progress information
)

# Load genome structure and sequence
genome = loader.load(
    gtf_file=Path("path/to/annotations.gtf"),
    fasta_file=Path("path/to/genome.fa.gz")
)

# Access genome components
print(f"Genome: {genome.name} ({genome.species})")
print(f"Chromosomes: {len(genome.chromosomes)}")
print(f"Genes: {len(genome.genes)}")

# Get a specific chromosome
chr1 = genome.chromosomes.get("chr1")
if chr1:
    print(f"Chromosome: {chr1.name}, Length: {chr1.length}")
    print(f"Genes on chr1: {len(chr1.genes)}")
    
    # Get sequence for a region
    region_seq = chr1.get_sequence(1000000, 1000100)
    print(f"Sequence: {region_seq}")

# Get a specific gene
tp53 = genome.genes.get("ENSG00000141510")
if tp53:
    print(f"TP53: {tp53.chrom}:{tp53.start}-{tp53.end} ({tp53.strand})")
    
    # Get gene sequence
    gene_seq = tp53.get_sequence()
    
    # Get coding sequence
    for transcript in tp53.transcripts:
        cds_seq = transcript.get_coding_sequence()
        protein = transcript.get_protein()
        print(f"Transcript {transcript.id}: CDS length: {len(cds_seq)}, Protein length: {len(protein)}")
```

### Working with MSI Sites

```python
from pathlib import Path
from pygnome.parsers.msi.msi_sites_reader import MsiSitesReader
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType

# Parse MSI sites file
reader = MsiSitesReader(Path("path/to/msi_sites.txt"))
msi_sites = reader.read_all()

# Create a specialized MSI store
msi_store = GenomicFeatureStore(store_type=StoreType.MSI)

# Add MSI sites to the store
with msi_store:
    for site in msi_sites:
        msi_store.add(site)

# Query MSI sites
sites_in_region = msi_store.get_by_interval("chr1", 1000000, 2000000)
for site in sites_in_region:
    print(f"MSI site: {site.chrom}:{site.start}-{site.end}, Repeat: {site.repeat_unit}")
```

## Performance Considerations

PyGnome offers multiple feature store implementations with different performance characteristics:

- **IntervalTreeStore**: Best for random access queries (default)
- **BinnedGenomicStore**: Good balance between memory usage and query speed
- **BruteForceFeatureStore**: Lowest memory usage but slower queries
- **MsiChromosomeStore**: Specialized for MSI sites

For large genomes, consider:

1. Using the context manager pattern when adding features to ensure proper indexing
2. Saving the populated store to disk with `store.save()` for faster loading in future sessions

Building genomic feature stores with large datasets can be time-consuming, especially when creating indexes for efficient querying. However, once built, these stores can be serialized to disk using Python's pickle format. This allows you to quickly load pre-built stores in future sessions, avoiding the need to rebuild them each time:

```python
# Save a populated store to disk (trimming is done automatically)
store.save(Path("path/to/store.pkl"))

# Later, quickly load the pre-built store
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.