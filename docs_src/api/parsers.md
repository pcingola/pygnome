# Parsers

The `parsers` module provides parsers for common genomic file formats. It includes parsers for FASTA/FASTQ, GFF/GTF, VCF, and MSI formats.

## Overview

The parsers module is designed to efficiently read and parse genomic data from various file formats. It provides a consistent interface for working with different file formats and converts the parsed data into PyGnome's object model.

Key features:

- Support for common genomic file formats (FASTA, FASTQ, GFF, GTF, VCF)
- Memory-efficient parsing of large files
- Conversion of parsed data into PyGnome's object model
- Support for compressed files (gzip, bgzip)

## FASTA/FASTQ Parsers

### FastaParser

```python
class FastaParser:
```

Parser for FASTA files. FASTA is a text-based format for representing nucleotide or peptide sequences.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the FASTA file

#### Methods

| Method | Description |
|--------|-------------|
| `load() -> list[FastaRecord]` | Load all records from the FASTA file |
| `load_as_dict() -> dict[str, FastaRecord]` | Load all records from the FASTA file as a dictionary |
| `parse_as_dna_strings(filepath: Path) -> dict[str, DnaString]` | Parse a FASTA file and return a dictionary of DnaString objects |

### FastqParser

```python
class FastqParser:
```

Parser for FASTQ files. FASTQ is a text-based format for storing both a nucleotide sequence and its corresponding quality scores.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the FASTQ file

#### Methods

| Method | Description |
|--------|-------------|
| `load() -> list[FastqRecord]` | Load all records from the FASTQ file |
| `load_as_dict() -> dict[str, FastqRecord]` | Load all records from the FASTQ file as a dictionary |

## GFF/GTF Parsers

### GffParser

```python
class GffParser:
```

Base class for GFF/GTF parsers. GFF (General Feature Format) and GTF (Gene Transfer Format) are tab-delimited text formats for describing genes and other features of DNA, RNA, and protein sequences.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the GFF/GTF file

#### Methods

| Method | Description |
|--------|-------------|
| `__iter__() -> Iterator[GffRecord]` | Iterate over records in the GFF/GTF file |
| `parse_attributes(attributes_str: str) -> dict[str, str]` | Parse the attributes field of a GFF/GTF record |

### Gff3Parser

```python
class Gff3Parser(GffParser):
```

Parser for GFF3 files. GFF3 is the latest version of the General Feature Format.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the GFF3 file

### GtfParser

```python
class GtfParser(GffParser):
```

Parser for GTF files. GTF (Gene Transfer Format) is a refinement of GFF that is used to describe genes and other features of DNA, RNA, and protein sequences.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the GTF file

## VCF Parsers

### VcfReader

```python
class VcfReader:
```

Parser for VCF files. VCF (Variant Call Format) is a text file format for storing gene sequence variations.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the VCF file

#### Methods

| Method | Description |
|--------|-------------|
| `__iter__() -> Iterator[VcfRecord]` | Iterate over records in the VCF file |
| `get_samples() -> list[str]` | Get the sample names from the VCF file |
| `get_header() -> VcfHeader` | Get the header from the VCF file |
| `fetch(chrom: str, start: int, end: int) -> Iterator[VcfRecord]` | Fetch records in a specific region |

### VcfRecord

```python
class VcfRecord:
```

Represents a single record (line) in a VCF file. Provides methods for accessing and modifying VCF fields.

#### Methods

| Method | Description |
|--------|-------------|
| `get_info(field_id: str) -> Any` | Get the value of an INFO field |
| `info[field_id]` | Dictionary-like access to INFO fields |
| `has_info(field_id: str) -> bool` | Check if an INFO field is present |
| `set_info(field_id: str, value: Any) -> None` | Set the value of an INFO field |
| `get_genotypes() -> list[Genotype]` | Get the genotypes for all samples |
| `get_genotype_value(field_id: str, sample_idx: int = 0) -> Any` | Get the value of a genotype field for a specific sample |
| `variant_annotations()` | Get a variant annotations parser for this record |
| `__iter__() -> Iterator[Variant]` | Iterate over the variants represented in this VCF record |

#### INFO Field Access

VcfRecord provides two ways to access INFO fields:

```python
# Method-based access
depth = record.get_info("DP")

# Dictionary-like access (more pythonic)
depth = record.info["DP"]

# Iterating through all INFO fields
for field_id, field_value in record.info:
    print(f"{field_id} = {field_value}")
```

#### Variant Annotations

VcfRecord provides a convenient method to access variant annotations:

```python
# Get variant annotations directly from the record
for vann in record.variant_annotations():
    print(f"Variant annotation: {vann}")
```

### VcfInfo

```python
class VcfInfo:
```

Class for handling INFO fields in VCF records. Provides methods for parsing, accessing, and modifying INFO fields in VCF records.

#### Methods

| Method | Description |
|--------|-------------|
| `get(field_id: str) -> Any` | Get the value of an INFO field |
| `__getitem__(field_id: str) -> Any` | Dictionary-like access to INFO fields |
| `has(field_id: str) -> bool` | Check if an INFO field is present |
| `set(field_id: str, value: Any) -> None` | Set the value of an INFO field |
| `remove(field_id: str) -> None` | Remove an INFO field |
| `field_ids() -> list[str]` | Get the names of all INFO fields |
| `__iter__()` | Iterate over the INFO field names and values |

#### Usage Examples

```python
# Access INFO fields using dictionary-like syntax
depth = vcf_record.info["DP"]

# Iterate through all INFO fields
for field_id, field_value in vcf_record.info:
    print(f"{field_id} = {field_value}")

# Set an INFO field
vcf_record.info["DP"] = 30
```
### AnnParser

```python
class AnnParser:
```

Parser for the ANN field in VCF records. The ANN field contains variant annotation information according to the VCF annotation format specification.

#### Constructor

```python
def __init__(self, record: VcfRecord):
```

- `record`: A VcfRecord object containing the ANN field

#### Methods

Method | Description |
|--------|-------------|
`__iter__() -> Iterator[VariantAnnotation]` | Iterate over annotations in the ANN field |
`parse() -> None` | Parse the ANN field from the record |

### EffectType

```python
class EffectType(str, Enum):
```

Enumeration of effect types for variant annotations. These types represent the specific effects a variant can have on genomic features. Each effect type has an associated impact level and can be mapped to a Sequence Ontology term.

#### Methods

| Method | Description |
|--------|-------------|
| `get_impact() -> AnnotationImpact` | Get the impact level for the effect type |
| `to_sequence_ontology() -> str` | Convert the effect type to a Sequence Ontology term |
| `from_sequence_ontology(so_term: str) -> EffectType | None` | Convert a Sequence Ontology term to an EffectType |

#### Impact Levels

Effect types are categorized into four impact levels:

- **HIGH**: Disruptive variants with high impact (e.g., frameshift, stop gained)
- **MODERATE**: Variants that might change protein effectiveness (e.g., missense)
- **LOW**: Variants unlikely to change protein behavior (e.g., synonymous)
- **MODIFIER**: Non-coding variants or variants affecting non-coding genes (e.g., intron)

#### Examples

```python
# Get the impact level of an effect type
effect = EffectType.FRAME_SHIFT
impact = effect.get_impact()  # Returns AnnotationImpact.HIGH

# Convert between effect types and Sequence Ontology terms
so_term = EffectType.NON_SYNONYMOUS_CODING.to_sequence_ontology()  # Returns "missense_variant"
effect = EffectType.from_sequence_ontology("missense_variant")  # Returns EffectType.NON_SYNONYMOUS_CODING
```

### VariantAnnotation

```python
class VariantAnnotation:
```

Represents a single annotation entry from the ANN field in a VCF record.

#### Constructor

```python
def __init__(self, allele: str, annotation: str, effect: EffectType, putative_impact: AnnotationImpact):
```

- `allele`: The allele being annotated
- `annotation`: The original annotation string from the VCF (e.g., "missense_variant")
- `effect`: The parsed effect type (e.g., EffectType.NON_SYNONYMOUS_CODING)
- `putative_impact`: The putative impact of the variant (HIGH, MODERATE, LOW, MODIFIER)

#### Properties

Property | Type | Description |
|----------|------|-------------|
`allele` | str | The allele being annotated |
`annotation` | str | The original annotation string from the VCF |
`effect` | EffectType | The parsed effect type |
`putative_impact` | AnnotationImpact | The putative impact of the variant |
`gene_name` | str | The gene name |
`gene_id` | str | The gene ID |
`feature_type` | FeatureType | The type of feature (e.g., transcript, motif) |
`feature_id` | str | The ID of the feature |
`transcript_biotype` | BiotypeCoding | The biotype of the transcript (Coding, Noncoding) |
`rank` | int | The rank of the exon or intron |
`total` | int | The total number of exons or introns |
`hgvs_c` | str | The HGVS notation at the DNA level |
`hgvs_p` | str | The HGVS notation at the protein level |
`cdna_pos` | int | The position in the cDNA |
`cdna_length` | int | The length of the cDNA |
`cds_pos` | int | The position in the CDS |
`cds_length` | int | The length of the CDS |
`protein_pos` | int | The position in the protein |
`protein_length` | int | The length of the protein |
`distance` | int | The distance to the feature |
`messages` | list[ErrorWarningType] | Error, warning, or information messages |

## MSI Parsers

### MsiSitesReader

```python
class MsiSitesReader:
```

Parser for MSI (Microsatellite Instability) sites files. MSI sites are regions of the genome with repetitive DNA sequences.

#### Constructor

```python
def __init__(self, filepath: Path):
```

- `filepath`: Path to the MSI sites file

#### Methods

| Method | Description |
|--------|-------------|
| `read_all() -> list[MsiSiteRecord]` | Read all MSI sites from the file |
| `read_by_chromosome(chrom: str) -> list[MsiSiteRecord]` | Read MSI sites for a specific chromosome |

## Genome Loader

### GenomeLoader

```python
class GenomeLoader:
```

Class for loading complete genomes from annotation and sequence files. This class combines sequence data from FASTA files with annotation data from GFF/GTF files to build a complete Genome object with chromosomes, genes, transcripts, exons, and other genomic features.

#### Constructor

```python
def __init__(self, annotation_file: Path, sequence_file: Path, genome_name: str = "genome", species: str = None, verbose: bool = False, error_handling: ErrorHandling = ErrorHandling.WARN):
```

- `annotation_file`: Path to the GFF/GTF annotation file
- `sequence_file`: Path to the FASTA sequence file
- `genome_name`: Name of the genome
- `species`: Species name
- `verbose`: Whether to print progress information during loading
- `error_handling`: How to handle consistency errors (throw, warn, or ignore)

#### Methods

| Method | Description |
|--------|-------------|
| `load() -> Genome` | Load a genome from the annotation and sequence files |
| `load_sequences(sequence_file: Path) -> dict[str, Chromosome]` | Load chromosome sequences from a FASTA file |
| `load_features(gff_file: Path) -> None` | Load genomic features from a GFF/GTF file |
| `check_consistency() -> None` | Check for consistency errors in the loaded genome |

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
        
        # Access INFO fields directly using dictionary-like syntax
        if record.has_info("DP"):
            depth = record.info["DP"]
            print(f"Read depth: {depth}")
        
        # Iterate through INFO fields
        for field_id, field_value in record.info:
            print(f"{field_id} = {field_value}")
        
        # Create variant objects from the record
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

# Open a VCF file
with VcfReader(Path("path/to/variants.vcf")) as reader:
    # Iterate through records
    for record in reader:
        # Get variant annotations directly from the record
        for vann in record.variant_annotations():
            print(f"Variant annotation: {vann.allele} - {vann.annotation}")
            print(f"  Impact: {vann.putative_impact}")
            
            if vann.gene_name:
                print(f"  Gene: {vann.gene_name}")
                
            if vann.feature_type and vann.feature_id:
                print(f"  Feature: {vann.feature_type.value} {vann.feature_id}")
                
            if vann.hgvs_c:
                print(f"  HGVS.c: {vann.hgvs_c}")
                
            if vann.hgvs_p:
                print(f"  HGVS.p: {vann.hgvs_p}")
```

### Loading a Complete Genome

```python
from pathlib import Path
from pygnome.parsers.genome_loader import GenomeLoader

# Create a genome loader
loader = GenomeLoader(
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz"),
    genome_name="GRCh38",
    species="Homo sapiens",
    verbose=True  # Print progress information
)

# Load genome structure and sequence
genome = loader.load()

# Access genome components
print(f"Genome: {genome.name} ({genome.species})")
print(f"Chromosomes: {len(genome.chromosomes)}")
print(f"Genes: {len(genome.genes)}")

# Get a specific chromosome
chr1 = genome.chromosomes.get("chr1")
if chr1:
    print(f"Chromosome: {chr1.name}, Length: {chr1.length}")
    print(f"Genes on chr1: {len(chr1.genes)}")
```