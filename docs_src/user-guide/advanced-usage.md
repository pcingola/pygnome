# Advanced Usage

This guide covers advanced usage patterns and optimization techniques for PyGnome.

## Performance Optimization

### Feature Store Selection

PyGnome offers multiple feature store implementations with different performance characteristics:

| Store Type | Memory Usage | Query Speed | Best For |
|------------|--------------|-------------|----------|
| `IntervalTreeStore` | Medium | Fast | General purpose, balanced performance |
| `BinnedGenomicStore` | Low | Medium | Large genomes, memory-constrained environments |
| `BruteForceFeatureStore` | Very Low | Slow | Testing, very small datasets |
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

### Memory Optimization

When working with large genomes, use these techniques to reduce memory usage:

#### 1. Use the Context Manager Pattern

Always use the context manager pattern when adding features to ensure proper indexing:

```python
with feature_store:
    for feature in features:
        feature_store.add(feature)
```

#### 2. Trim Feature Stores Before Serialization

Call `trim()` to reduce memory usage before saving:

```python
# Trim the store to reduce memory usage
feature_store.trim()

# Save the trimmed store
feature_store.save(Path("path/to/store.pkl"))
```

#### 3. Use DnaStringArray for Multiple Sequences

When working with many small sequences, use `DnaStringArray` instead of multiple `DnaString` objects:

```python
from pygnome.sequences.dna_string_array import DnaStringArray

# Create a DNA string array
array = DnaStringArray()

# Add sequences
for seq in sequences:
    array.add(seq)

# Trim to reduce memory usage
array.trim()
```

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

### Parallel Processing

Use Python's multiprocessing to process chromosomes in parallel:

```python
import multiprocessing as mp
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader

def process_chromosome(chrom_name, genome):
    # Get the chromosome
    chrom = genome.chromosomes.get(chrom_name)
    if not chrom:
        return None
    
    # Process the chromosome
    # ...
    
    return result

# Load the genome
loader = GenomeLoader(genome_name="GRCh38", species="Homo sapiens")
genome = loader.load(
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz")
)

# Process chromosomes in parallel
with mp.Pool(processes=mp.cpu_count()) as pool:
    results = pool.starmap(
        process_chromosome,
        [(chrom_name, genome) for chrom_name in genome.chromosomes.keys()]
    )
```

## Advanced Feature Store Usage

### Custom Feature Filtering

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

## Working with Genomic Variants

### Annotating Variants with Genomic Features

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
        # Get the variant position
        chrom = record.get_chrom()
        pos = record.get_pos()
        
        # Find features at this position
        features = store.get_by_position(chrom, pos)
        
        # Annotate the variant
        for feature in features:
            print(f"Variant at {chrom}:{pos} overlaps {feature.id} ({feature.__class__.__name__})")
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