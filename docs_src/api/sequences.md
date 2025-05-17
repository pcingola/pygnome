# Sequences

The `sequences` module provides memory-efficient representations of DNA and RNA sequences. It includes classes for representing individual sequences and collections of sequences.

## Overview

The sequences module is designed to efficiently store and manipulate nucleotide sequences. It uses a 2-bit encoding for DNA and RNA sequences, which allows for significant memory savings compared to string-based representations.

Key features:

- Memory-efficient 2-bit representation (A=00, C=01, G=10, T/U=11)
- Support for common sequence operations (complement, reverse complement, transcription, translation)
- Efficient storage of multiple sequences in a single array
- Support for slicing and indexing

## Core Classes

### BaseSequence

```python
class BaseSequence(ABC):
```

Abstract base class for efficient 2-bit representation of nucleotide sequences. This class provides common functionality for storing nucleotide sequences using 2 bits per nucleotide, allowing 16 nucleotides to be packed into a single 32-bit integer.

#### Constructor

```python
def __init__(self, sequence: str):
```

- `sequence`: A string containing nucleotides

#### Methods

| Method | Description |
|--------|-------------|
| `__len__() -> int` | Return the length of the sequence |
| `__str__() -> str` | Return the sequence as a string |
| `__getitem__(key) -> str` | Get a nucleotide or subsequence |
| `to_string() -> str` | Convert the entire sequence to a string |
| `substring(start: int, length: int \| None = None) -> str` | Extract a substring from the sequence |
| `__eq__(other) -> bool` | Check if two sequence objects are equal |

### DnaString

```python
class DnaString(BaseSequence):
```

Efficient 2-bit representation of DNA sequences. This class stores DNA sequences (A, C, G, T) using 2 bits per nucleotide, allowing 16 nucleotides to be packed into a single 32-bit integer.

#### Constructor

```python
def __init__(self, sequence: str):
```

- `sequence`: A string containing DNA nucleotides (A, C, G, T)

#### Methods

| Method | Description |
|--------|-------------|
| `complement() -> DnaString` | Return the complement of the sequence |
| `reverse_complement() -> DnaString` | Return the reverse complement of the sequence |
| `transcribe() -> RnaString` | Transcribe DNA to RNA |
| `gc_content() -> float` | Calculate the GC content of the sequence |

### RnaString

```python
class RnaString(BaseSequence):
```

Efficient 2-bit representation of RNA sequences. This class stores RNA sequences (A, C, G, U) using 2 bits per nucleotide, allowing 16 nucleotides to be packed into a single 32-bit integer.

#### Constructor

```python
def __init__(self, sequence: str):
```

- `sequence`: A string containing RNA nucleotides (A, C, G, U)

#### Methods

| Method | Description |
|--------|-------------|
| `complement() -> RnaString` | Return the complement of the sequence |
| `reverse_complement() -> RnaString` | Return the reverse complement of the sequence |
| `translate() -> str` | Translate RNA to protein |

### DnaStringArray

```python
class DnaStringArray:
```

Efficient storage for millions of small DNA strings in a single NumPy array. This class stores multiple DNA sequences using the same 2-bit encoding as DnaString, but packs all sequences into a single contiguous array for improved memory efficiency when dealing with large numbers of sequences.

#### Constructor

```python
def __init__(self, initial_data_bytes: int = DEFAULT_CAPACITY, initial_strings: int = DEFAULT_NUMBER_OF_STRINGS):
```

- `initial_data_bytes`: Initial capacity in bytes for the data array (default: 1MB)
- `initial_strings`: Initial number of strings the array can hold (default: 100,000)

#### Methods

| Method | Description |
|--------|-------------|
| `add(sequence: str) -> int` | Add a DNA sequence to the array |
| `add_multiple(sequences: list[str]) -> list[int]` | Add multiple DNA sequences to the array |
| `get(idx: int) -> str` | Get a sequence by its index |
| `get_subsequence(idx: int, start: int, length: int \| None = None) -> str` | Extract a subsequence from a sequence in the array |
| `get_length(idx: int) -> int` | Get the length of a sequence |
| `__getitem__(idx: int) -> str` | Get a sequence by its index |
| `__len__() -> int` | Return the number of sequences in the array |
| `trim() -> None` | Trim the internal arrays to their actual used size |
| `get_stats() -> tuple[int, int, float]` | Get statistics about memory usage |
| `to_dna_string(idx: int) -> DnaString` | Convert a sequence in the array to a DnaString object |

## Usage Examples

### Working with DNA Sequences

```python
from pygnome.sequences.dna_string import DnaString

# Create a DNA sequence
dna = DnaString("ATGCATGCATGC")
print(f"Length: {len(dna)}")

# Get a subsequence
subseq = dna[3:9]  # Returns a new DnaString
print(f"Subsequence: {subseq}")

# Complement and reverse complement
comp = dna.complement()
rev_comp = dna.reverse_complement()
print(f"Complement: {comp}")
print(f"Reverse complement: {rev_comp}")

# Transcribe DNA to RNA
rna = dna.transcribe()  # Returns an RnaString
print(f"RNA: {rna}")
```

### Working with RNA Sequences

```python
from pygnome.sequences.rna_string import RnaString

# Create an RNA sequence
rna = RnaString("AUGCAUGCAUGC")
print(f"Length: {len(rna)}")

# Get a subsequence
subseq = rna[3:9]  # Returns a new RnaString
print(f"Subsequence: {subseq}")

# Complement and reverse complement
comp = rna.complement()
rev_comp = rna.reverse_complement()
print(f"Complement: {comp}")
print(f"Reverse complement: {rev_comp}")

# Translate RNA to protein
protein = rna.translate()
print(f"Protein: {protein}")
```

### Working with Multiple Sequences

```python
from pygnome.sequences.dna_string_array import DnaStringArray

# Create a DNA string array
array = DnaStringArray()

# Add sequences
idx1 = array.add("ATGCATGC")
idx2 = array.add("GCTAGCTA")

# Add multiple sequences at once
indices = array.add_multiple(["AAAAAA", "CCCCCC", "GGGGGG"])

# Access sequences
seq1 = array[idx1]
print(f"Sequence 1: {seq1}")

# Get subsequences
subseq = array.get_subsequence(idx2, 2, 4)
print(f"Subsequence: {subseq}")

# Get statistics
count, total_nt, bits_per_nt = array.get_stats()
print(f"Sequences: {count}, Total nucleotides: {total_nt}, Bits per nucleotide: {bits_per_nt:.2f}")

# Trim to reduce memory usage
array.trim()
```

## Memory Efficiency

The 2-bit encoding used by the sequences module provides significant memory savings compared to string-based representations:

- A standard Python string uses 1 byte (8 bits) per character
- DnaString and RnaString use 2 bits per nucleotide
- This results in a 4x reduction in memory usage

For example, a 1 million base pair chromosome would use:

- String representation: ~1 MB
- DnaString representation: ~250 KB

The DnaStringArray class provides even greater memory efficiency when storing multiple sequences by:

1. Minimizing memory overhead compared to individual DnaString objects
2. Improving cache locality for faster access patterns
3. Reducing memory fragmentation

## Performance Considerations

- Use DnaStringArray instead of multiple DnaString objects when working with many small sequences
- Call trim() on DnaStringArray before serialization to reduce memory usage
- For very large sequences (e.g., entire chromosomes), consider using memory-mapped files or chunked processing