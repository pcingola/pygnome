# Basic Usage

This guide covers the basic usage of PyGnome's core components with practical examples.

## Working with Genomic Models

PyGnome provides a comprehensive object model for genomic features.

### Creating a Genome

```python
from pygnome.genomics.genome import Genome
from pygnome.genomics.chromosome import Chromosome
from pygnome.genomics.gene import Gene
from pygnome.genomics.transcript import Transcript
from pygnome.genomics.exon import Exon
from pygnome.genomics.strand import Strand

# Create a genome
genome = Genome(name="MyGenome", species="Example species")

# Create a chromosome
chromosome = Chromosome(name="chr1", length=248956422)

# Add the chromosome to the genome
genome.add_chromosome(chromosome)

# Create a gene
gene = Gene(
    id="GENE001",
    name="Example Gene",
    chrom="chr1",
    start=1000,
    end=5000,
    strand=Strand.POSITIVE
)

# Create a transcript
transcript = Transcript(
    id="TRANSCRIPT001",
    gene_id="GENE001",
    chrom="chr1",
    start=1000,
    end=5000,
    strand=Strand.POSITIVE
)

# Create exons
exon1 = Exon(
    id="EXON001",
    chrom="chr1",
    start=1000,
    end=2000,
    strand=Strand.POSITIVE
)

exon2 = Exon(
    id="EXON002",
    chrom="chr1",
    start=3000,
    end=5000,
    strand=Strand.POSITIVE
)

# Add exons to transcript
transcript.exons = [exon1, exon2]

# Add transcript to gene
gene.add_transcript(transcript)

# Add gene to chromosome
chromosome.add_gene(gene)

# Now the gene is also accessible from the genome
print(genome.genes["GENE001"])
```

## Parsing Genomic Files

PyGnome provides parsers for common genomic file formats.

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
        
        # Create variant objects from the record
        for variant in record:  # Uses VariantFactory internally
            print(f"Variant: {variant}")
            
        # Access genotypes
        genotypes = record.get_genotypes()
        for i, genotype in enumerate(genotypes):
            print(f"  {samples[i]}: {genotype}")
```

## Working with DNA/RNA Sequences

PyGnome provides memory-efficient representations of DNA and RNA sequences.

### Creating and Manipulating DNA Sequences

```python
from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString

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

# Create an RNA sequence
rna = RnaString("AUGCAUGCAUGC")

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
```

## Using Feature Stores

Feature stores provide efficient storage and retrieval of genomic features.

### Creating and Using a Feature Store

```python
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.genomics.gene import Gene
from pygnome.genomics.strand import Strand
from pathlib import Path

# Create a feature store using interval trees (default)
store = GenomicFeatureStore()

# Or choose a different implementation
binned_store = GenomicFeatureStore(store_type=StoreType.BINNED, bin_size=100000)
brute_force_store = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)

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

# Save and load the store
store.save(Path("path/to/store.pkl"))
loaded_store = GenomicFeatureStore.load(Path("path/to/store.pkl"))
```

## Loading a Complete Genome

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
    annotation_file=Path("path/to/annotations.gtf"),
    sequence_file=Path("path/to/genome.fa.gz")
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

# Get a specific gene
tp53 = genome.genes.get("TP53")
if tp53:
    print(f"TP53: {tp53.chrom}:{tp53.start}-{tp53.end} ({tp53.strand})")
    
    # Get gene transcripts
    for transcript in tp53.transcripts:
        print(f"Transcript {transcript.id}: Exons: {len(transcript.exons)}")