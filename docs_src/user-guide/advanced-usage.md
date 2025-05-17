# Advanced Usage

This guide covers advanced usage patterns and optimization techniques for PyGnome.

## Genomic Feature Stores

### Genomic Feature Store Selection

Genomic feature stores are one of the core solutions in PyGnome, providing specialized data structures for efficient storage, indexing, and querying of genomic features based on their genomic coordinates. They solve the fundamental bioinformatics challenge of quickly locating genomic elements within large genomes.

PyGnome offers multiple feature store implementations with different performance characteristics:

| Store Type | Memory Usage | Query Speed | Best For |
|------------|--------------|-------------|----------|
| `IntervalTreeStore` | Medium | Fast | General purpose, balanced performance |
| `BinnedGenomicStore` | Medium | Medium | Large genomes, memory-constrained environments |
| `BruteForceFeatureStore` | Medium | Slow | Testing, very small datasets |
| `MsiChromosomeStore` | Very Low | Fast | Specialized for microsatellite sites |

Choose the appropriate store type based on your specific use case:

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType

# For general use
default_store = GenomicFeatureStore()  # Uses IntervalTreeStore

# For memory-constrained environments
binned_store = GenomicFeatureStore(store_type=StoreType.BINNED, bin_size=100000)

# For testing with small datasets
brute_force_store = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)

# For microsatellite sites
msi_store = GenomicFeatureStore(store_type=StoreType.MSI)
```

### Memory and Performance Optimization

When working with large genomes, use these techniques to optimize performance:

#### 1. Use the Context Manager Pattern

Always use the context manager pattern when adding features to ensure proper indexing:

```python
with feature_store:
    for feature in features:
        feature_store.add(feature)
```

#### 2. Save and Load Feature Stores to Avoid Rebuilding

Building genomic feature stores with large datasets can be time-consuming, especially when creating indexes for efficient querying. Once built, save them to disk to avoid rebuilding them in future sessions:

```python
# Building a store with millions of features can take time
store = GenomicFeatureStore()
with store:
    # Adding features from a large genome...
    for gene in genome.genes.values():
        store.add(gene)
        # Add transcripts, exons, etc.

# Save the built store (trimming is done automatically)
store.save(Path("path/to/store.pkl"))

# In future sessions, quickly load the pre-built store
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
# Ready to use immediately without rebuilding indexes
```

Note: The `save()` method automatically calls `trim()` to reduce memory usage before serialization, so you don't need to call it manually.

## Working with Large Genomes

When working with large genomes, consider these approaches:

### Chunked Processing

Process large files in chunks to reduce memory usage:

```python
from pathlib import Path
from pygnome.parsers.fasta.fasta_parser import FastaParser

# Process a large FASTA file in chunks
parser = FastaParser(Path("path/to/large_genome.fa"))
for record in parser:
    # Process one sequence at a time
    process_sequence(record.identifier, record.sequence)
```

## Advanced Feature Store Usage

## Working with Genomic Variants

### Reading Variants from VCF Files

```python
from pathlib import Path
from pygnome.parsers.vcf.vcf_reader import VcfReader

# Open a VCF file
with VcfReader(Path("path/to/variants.vcf")) as reader:
    for vcf_record in reader:
        # A single VCF record (i.e. vcf line) can have multiple variants
        for variant in vcf_record:
            print(f"Variant: {variant}")
```

### Annotating Variants with Genomic Feature stores

```python
from pathlib import Path
from pygnome.parsers.vcf.vcf_reader import VcfReader
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore
from pygnome.parsers.genome_loader import GenomeLoader

# Load the genome
loader = GenomeLoader(genome_name="GRCh38", species="Homo sapiens")
genome = loader.load(
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz")
)

# Create a feature store with all genes
store = GenomicFeatureStore()
with store:
    for gene in genome.genes.values():
        store.add(gene)
        for transcript in gene.transcripts:
            store.add(transcript)
            for exon in transcript.exons:
                store.add(exon)

# Open a VCF file
with VcfReader(Path("path/to/variants.vcf")) as reader:
    for record in reader:
        # A single VCF record (i.e. vcf line) can have multiple variants
        for variant in vcf_record:
            # Get the variant position
            chrom = variant.get_chrom()
            pos = variant.get_pos()

            # Find features at this position
            features = store.get_by_position(chrom, pos)
            
            # Show matches
            for feature in features:
                print(f"Variant at {chrom}:{pos} overlaps {feature.id} ({feature.__class__.__name__})")
```

## Custom Feature Filtering

Combine feature store queries with custom filtering:

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore
from pygnome.genomics.gene import Gene

# Create and populate a feature store
store = GenomicFeatureStore()
# ... add features ...

# Get all features in a range
features = store.get_by_interval("chr1", 1000000, 2000000)

# Filter for specific feature types
genes = [f for f in features if isinstance(f, Gene)]

# Filter by additional criteria
protein_coding_genes = [g for g in genes if g.biotype == "protein_coding"]
```

## Custom Genomic Feature Types

You can create custom genomic feature types by extending the `GenomicFeature` class:

```python
from dataclasses import dataclass
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand

@dataclass
class EnhancerRegion(GenomicFeature):
    """A genomic enhancer region."""
    target_genes: list[str] = None
    activity_score: float = 0.0
    
    def __post_init__(self):
        super().__post_init__()
        if self.target_genes is None:
            self.target_genes = []
    
    def is_active(self, threshold: float = 0.5) -> bool:
        """Check if the enhancer is active based on its activity score."""
        return self.activity_score >= threshold

# Create an enhancer
enhancer = EnhancerRegion(
    id="ENH001",
    chrom="chr1",
    start=1000000,
    end=1001000,
    strand=Strand.POSITIVE,
    target_genes=["GENE001", "GENE002"],
    activity_score=0.75
)

# Add to a feature store
store = GenomicFeatureStore()
with store:
    store.add(enhancer)

# Query enhancers
enhancers = [f for f in store.get_by_interval("chr1", 900000, 1100000) 
             if isinstance(f, EnhancerRegion)]
active_enhancers = [e for e in enhancers if e.is_active()]
```