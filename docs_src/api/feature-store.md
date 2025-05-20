# Genomic Feature Store

The `feature_store` module provides efficient storage and retrieval of genomic features. It offers multiple implementations with different performance characteristics to suit various use cases.

## Overview

*Genomic feature stores* are one of the core solutions in PyGnome, providing specialized data structures for efficient storage, indexing, and querying of genomic features based on their genomic coordinates. They solve the fundamental bioinformatics challenge of quickly locating genomic elements within large genomes.

### What is a Genomic Feature Store?

A genomic feature store is a data structure that:

1. **Stores genomic features** (genes, transcripts, exons, variants, etc.) organized by their chromosomal locations
2. **Indexes these features** using efficient spatial data structures (interval trees, binning, etc.)
3. **Provides fast query operations** to find features based on genomic coordinates
4. **Optimizes memory usage** through various implementation strategies

The feature store is designed to efficiently store and query genomic features based on their genomic coordinates. It supports several types of queries:

- **Position queries**: Find all features at a specific position
- **Interval queries**: Find all features that overlap with a given range
- **Nearest feature queries**: Find the nearest feature to a specific position
- **Variant queries**: Find features that match a specific genomic variant (chr, pos, ref, alt)

## Core Components

### GenomicFeatureStoreProtocol

```python
class GenomicFeatureStoreProtocol(Protocol):
```

This protocol defines the interface that all feature store implementations must follow. It includes methods for adding features, querying features, and managing the store.

#### Methods

| Method | Description |
|--------|-------------|
| `add(feature: GenomicFeature) -> None` | Add a genomic feature to the store |
| `get_by_position(chrom: str, position: int) -> list[GenomicFeature]` | Get all features at a specific position |
| `get_by_interval(chrom: str, start: int, end: int) -> list[GenomicFeature]` | Get all features that overlap with the given range |
| `get_nearest(chrom: str, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature \| None` | Get the nearest feature to the given position |
| `get_by_variant(chrom: str, position: int, ref: str, alt: str) -> list[GenomicFeature]` | Get all features that match the specific variant (chr, pos, ref, alt) |
| `get_variant(variant: Variant) -> list[GenomicFeature]` | Get all features that match the given Variant object |
| `__getitem__(chrom: str) -> ChromosomeFeatureStore` | Get a chromosome store by name |
| `__iterator__()` | Iterate over all chromosome stores |
| `trim() -> None` | Trim internal data structures to reduce memory usage |

### GenomicFeatureStore

```python
class GenomicFeatureStore(GenomicFeatureStoreProtocol):
```

The main implementation of the genomic feature store. It delegates storage and queries to chromosome-specific stores based on the chosen store type.

#### Constructor

```python
def __init__(self, store_type: StoreType | str = StoreType.INTERVAL_TREE, bin_size: int = 100000):
```

- `store_type`: Type of store to use (default: `StoreType.INTERVAL_TREE`)
- `bin_size`: Size of bins for the binned store (default: 100000)

#### Methods

| Method | Description |
|--------|-------------|
| `add(feature: GenomicFeature) -> None` | Add a genomic feature to the store |
| `add_features(features: list[GenomicFeature]) -> None` | Add multiple genomic features to the store |
| `get_by_position(chrom: str, position: int) -> list[GenomicFeature]` | Get all features at a specific position |
| `get_by_interval(chrom: str, start: int, end: int) -> list[GenomicFeature]` | Get all features that overlap with the given range |
| `get_nearest(chrom: str, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature \| None` | Get the nearest feature to the given position |
| `get_by_variant(chrom: str, position: int, ref: str, alt: str) -> list[GenomicFeature]` | Get all features that match the specific variant (chr, pos, ref, alt) |
| `get_variant(variant: Variant) -> list[GenomicFeature]` | Get all features that match the given Variant object |
| `get_chromosomes() -> list[str]` | Get all chromosome names in the store |
| `get_feature_count(chrom: str \| None = None) -> int` | Get the number of features in the store |
| `trim() -> None` | Trim internal data structures to reduce memory usage |
| `save(filepath: Path) -> None` | Save the genomic feature store to a file using pickle |
| `load(filepath: Path) -> GenomicFeatureStore` | Load a genomic feature store from a file |

#### Context Manager

The `GenomicFeatureStore` class implements the context manager protocol (`__enter__` and `__exit__` methods), which should be used when adding features to ensure proper indexing:

```python
with store:
    for feature in features:
        store.add(feature)
```

### StoreType

```python
class StoreType(str, Enum):
```

An enumeration of the available feature store types.

| Value | Description |
|-------|-------------|
| `INTERVAL_TREE` | Uses interval trees for efficient range queries |
| `BINNED` | Uses binning for memory-efficient storage |
| `BRUTE_FORCE` | Simple implementation for testing |
| `MSI` | Specialized for microsatellite instability sites |

### ChromosomeFeatureStore

```python
class ChromosomeFeatureStore(ABC):
```

Abstract base class for chromosome-specific genomic feature storage. It has a list of features and provides methods to add and query them.

#### Constructor

```python
def __init__(self, chromosome: str) -> None:
```

- `chromosome`: Name of the chromosome

#### Methods

| Method | Description |
|--------|-------------|
| `add(feature: GenomicFeature) -> None` | Add a feature to this chromosome's store |
| `get_by_position(position: int) -> list[GenomicFeature]` | Get all features at a specific position |
| `get_by_interval(start: int, end: int) -> list[GenomicFeature]` | Get all features that overlap with the given range |
| `get_features() -> list[GenomicFeature]` | Get all features |
| `get_nearest(position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature \| None` | Get the nearest feature to the given position |
| `get_by_variant(position: int, ref: str, alt: str) -> list[GenomicFeature]` | Get all features that match the specific variant (pos, ref, alt) |
| `get_variant(variant: Variant) -> list[GenomicFeature]` | Get all features that match the given Variant object |
| `index_build_start() -> None` | Start building the index |
| `index_build_end() -> None` | Finish building the index |
| `trim() -> None` | Trim internal data structures to reduce memory usage |

## Store Implementations

### IntervalTreeStore

```python
class IntervalTreeStore(ChromosomeFeatureStore):
```

Store genomic features using an efficient interval tree. This is the default implementation and provides a good balance between memory usage and query speed.

#### Performance Characteristics

- Memory Usage: Medium
- Query Speed: Fast
- Best For: General purpose, balanced performance

### BinnedGenomicStore

```python
class BinnedGenomicStore(ChromosomeFeatureStore):
```

Store genomic features using a memory-efficient binning approach. Features are grouped into bins based on their genomic coordinates, which allows for efficient range queries while using less memory than interval trees.

#### Constructor

```python
def __init__(self, chromosome: str, bin_size: int = DEFAULT_BIN_SIZE):
```

- `chromosome`: Name of the chromosome
- `bin_size`: Size of each bin in base pairs (default: 100000)

#### Performance Characteristics

- Memory Usage: Low
- Query Speed: Medium
- Best For: Large genomes, memory-constrained environments

### BruteForceFeatureStore

```python
class BruteForceFeatureStore(ChromosomeFeatureStore):
```

A naive brute-force implementation for genomic feature storage. This is not memory efficient and is not recommended for large datasets. It is primarily for testing purposes.

#### Performance Characteristics

- Memory Usage: Very Low
- Query Speed: Slow
- Best For: Testing, very small datasets

### MsiChromosomeStore

```python
class MsiChromosomeStore(ChromosomeFeatureStore):
```

Efficient storage for millions of MSI (Microsatellite Instability) sites in a chromosome. Uses NumPy arrays and DnaStringArray for memory efficiency. Implements binary search for efficient querying.

#### Constructor

```python
def __init__(self, chrom: str, feature_count: int = DEFAULT_FEATURE_COUNT, max_lengths_by_bin: dict[int, int] | None = None, bin_size: int = DEFAULT_BIN_SIZE):
```

- `chrom`: Name of the chromosome
- `feature_count`: Number of features to allocate space for (default: 1024)
- `max_lengths_by_bin`: Dictionary mapping bin IDs to maximum feature length in that bin
- `bin_size`: Size of each bin in base pairs (default: 100000)

#### Performance Characteristics

- Memory Usage: Very Low
- Query Speed: Fast
- Best For: Specialized for microsatellite sites

## Usage Examples

### Basic Usage

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore
from pygnome.genomics.gene import Gene
from pygnome.genomics.strand import Strand

# Create a feature store
store = GenomicFeatureStore()

# Create some features
gene1 = Gene(id="GENE001", chrom="chr1", start=1000, end=5000, strand=Strand.POSITIVE)
gene2 = Gene(id="GENE002", chrom="chr1", start=7000, end=9000, strand=Strand.NEGATIVE)
gene3 = Gene(id="GENE003", chrom="chr2", start=2000, end=6000, strand=Strand.POSITIVE)

# Add features to the store
with store:  # Use context manager to ensure proper indexing
    store.add(gene1)
    store.add(gene2)
    store.add(gene3)

# Query features
features_at_position = store.get_by_position("chr1", 1500)
print(f"Features at position chr1:1500: {features_at_position}")

features_in_range = store.get_by_interval("chr1", 4000, 8000)
print(f"Features in range chr1:4000-8000: {features_in_range}")

nearest_feature = store.get_nearest("chr1", 6000)
print(f"Nearest feature to chr1:6000: {nearest_feature}")

# Query by variant
from pygnome.genomics.variant import Variant, SNP

# Query by variant parameters
variant_features = store.get_by_variant("chr1", 3000, "A", "G")
print(f"Features matching variant chr1:3000 A>G: {variant_features}")

# Or query using a Variant object
snp = SNP(id="rs123", chrom="chr1", start=3000, end=3001, strand=Strand.POSITIVE, ref="A", alt="G")
variant_features = store.get_variant(snp)
print(f"Features matching variant object: {variant_features}")
```

### Choosing a Store Type

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType

# Create a feature store with interval trees (default)
default_store = GenomicFeatureStore()

# Create a feature store with binning
binned_store = GenomicFeatureStore(store_type=StoreType.BINNED, bin_size=100000)

# Create a feature store with brute force
brute_force_store = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)

# Create a feature store for MSI sites
msi_store = GenomicFeatureStore(store_type=StoreType.MSI)
```

### Saving and Loading

```python
from pathlib import Path
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore

# Create and populate a feature store
store = GenomicFeatureStore()
# ... add features ...

# Trim the store to reduce memory usage
store.trim()

# Save the store
store.save(Path("path/to/store.pkl"))

# Load the store
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
```

## Performance Considerations

### Memory Usage

The memory usage of the feature store depends on the implementation:

- `IntervalTreeStore`: Uses more memory but provides faster queries
- `BinnedGenomicStore`: Uses less memory but queries may be slightly slower
- `BruteForceFeatureStore`: Uses minimal memory but queries are slow
- `MsiChromosomeStore`: Specialized for MSI sites, very memory efficient

### Query Speed

The query speed depends on the implementation and the number of features:

- `IntervalTreeStore`: O(log n + k) for finding k intervals that overlap with a given point or range
- `BinnedGenomicStore`: O(b + k) where b is the number of bins that overlap with the query range
- `BruteForceFeatureStore`: O(n) where n is the total number of features
- `MsiChromosomeStore`: O(log n + k) using binary search

### Best Practices

1. Use the context manager pattern when adding features to ensure proper indexing
2. Choose the appropriate store type based on your specific use case
3. For large genomes, consider saving the populated store to disk with `store.save()` for faster loading in future sessions

### Build Time vs. Load Time

Building genomic feature stores with large datasets can be time-consuming, especially when creating indexes for efficient querying. The build time depends on:

- The number of features being added
- The complexity of the indexing structure (interval trees, binning, etc.)
- The performance characteristics of the chosen store implementation

However, once built, these stores can be serialized to disk using Python's pickle format. This allows you to quickly load pre-built stores in future sessions, avoiding the need to rebuild them each time:

```python
# Building a store can be time-consuming
store = GenomicFeatureStore()
with store:
    # Adding thousands or millions of features...
    for gene in genome.genes.values():
        store.add(gene)
        for transcript in gene.transcripts:
            store.add(transcript)
            # ...and so on

# Save the built store to avoid rebuilding it next time
# Note: trimming is done automatically during save
store.save(Path("path/to/store.pkl"))

# In future sessions, quickly load the pre-built store
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
# Ready to use immediately without rebuilding indexes
```

This approach is particularly valuable in production environments or when working with large reference genomes where the feature set doesn't change frequently.