# Py-gnome: Python library for genome annotations

A Python library for genomic annotations, similar to SnpEff/SnpSift but as a Python library.

## Features

### Efficient DNA/RNA Sequence Representation

- `DnaString`: Efficient 2-bit representation of DNA sequences (A, C, G, T)
- `RnaString`: Efficient 2-bit representation of RNA sequences (A, C, G, U)

Both classes provide:
- Memory-efficient storage (16 nucleotides per 32-bit integer)
- String-like interface (indexing, slicing, etc.)
- Substring extraction

### Genomic File Format Parsers

- `FastaParser`: Parser for FASTA format files
- `FastqParser`: Parser for FASTQ format files

Both parsers provide:
- Context manager interface for efficient file handling
- Iterator interface for streaming large files
- Utility methods for common operations (e.g., parsing as dictionary)
- Support for converting sequences to efficient `DnaString`/`RnaString` representations

## Usage Examples

### Working with DNA Sequences

```python
from pygnome.sequences import DnaString

# Create a DNA sequence
dna = DnaString("ACGTACGTACGTACGT")

# Access individual nucleotides
print(dna[0])  # "A"

# Slice the sequence
print(dna[1:5])  # "CGTA"

# Get a substring
print(dna.substring(2, 4))  # "GTAC"
```

### Working with RNA Sequences

```python
from pygnome.sequences import RnaString

# Create an RNA sequence (T's are automatically converted to U's)
rna = RnaString("ACGTACGTACGTACGT")

# Access individual nucleotides
print(rna[0])  # "A"
print(rna[3])  # "U" (not "T")
```

### Parsing FASTA Files

```python
from pathlib import Path
from pygnome.parsers import FastaParser

# Parse a FASTA file
fasta_file = Path("path/to/sequences.fasta")
for record in FastaParser.parse(fasta_file):
    print(f"Sequence {record.identifier}: {record.sequence[:10]}...")

# Parse as dictionary
sequences = FastaParser.parse_as_dict(fasta_file)
print(sequences["seq1"])

# Parse as DnaString objects
dna_sequences = FastaParser.parse_as_dna_strings(fasta_file)
print(dna_sequences["seq1"])
```

### Parsing FASTQ Files

```python
from pathlib import Path
from pygnome.parsers import FastqParser

# Parse a FASTQ file
fastq_file = Path("path/to/reads.fastq")
for record in FastqParser.parse(fastq_file):
    print(f"Read {record.identifier}: {record.sequence[:10]}...")
    print(f"Quality: {record.quality[:10]}...")

# Get Phred quality scores
for record in FastqParser.parse(fastq_file):
    quality_scores = record.get_quality_scores()
    print(f"Average quality: {sum(quality_scores) / len(quality_scores):.2f}")
```

## Development

### Running Tests

```bash
cd pygnome
python -m unittest discover -s tests