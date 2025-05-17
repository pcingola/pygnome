# VCF Reader Implementation Design

Based on the VCF format specification, this document outlines the implementation for a Python library to handle VCF files with lazy parsing techniques. The implementation will use zero-based indexing and half-open intervals as required.

## Core Classes

### 1. `VcfReader`
The main class responsible for reading and parsing VCF files.

- **Responsibilities**:
  - Open and read VCF files
  - Parse header information
  - Iterate through records
  - Provide access to metadata
  - Support indexing for random access

- **Key methods**:
  - `__init__(file_path)`: Initialize with file path
  - `parse_header()`: Parse header lines
  - `__iter__()`: Iterate through records
  - `fetch(chrom, start, end)`: Get records in a region (requires indexing)

### 2. `VcfHeader`
Represents the VCF file header.

- **Responsibilities**:
  - Store meta-information lines
  - Parse and provide access to metadata (contigs, INFO fields, FORMAT fields, etc.)
  - Validate record fields against header definitions

- **Key methods**:
  - `get_info_field_definition(id)`: Get definition for an INFO field
  - `get_format_field_definition(id)`: Get definition for a FORMAT field
  - `get_contigs()`: Get list of contigs

### 3. `VcfRecord`
Represents a single VCF record with lazy parsing.

- **Responsibilities**:
  - Store raw record data
  - Provide access to fixed fields (CHROM, POS, etc.)
  - Lazily parse INFO and FORMAT fields only when requested
  - Convert between 1-based VCF positions and 0-based internal positions

- **Key methods**:
  - `get_chrom()`, `get_pos()`, etc.: Access fixed fields
  - `get_info(field_id)`: Lazily parse and return specific INFO field
  - `get_format(field_id, sample_idx)`: Lazily parse and return FORMAT field for a sample
  - `get_genotypes()`: Lazily parse and return genotype data

### 4. `InfoField`
Represents an INFO field with type-specific parsing.

- **Responsibilities**:
  - Parse INFO field values according to their type
  - Convert between string representation and Python types

- **Key methods**:
  - `parse(value, field_type, number)`: Parse value according to type and number

### 5. `FormatField`
Represents a FORMAT field with type-specific parsing.

- **Responsibilities**:
  - Parse FORMAT field values according to their type
  - Handle per-sample data

- **Key methods**:
  - `parse(values, field_type, number)`: Parse values according to type and number

### 6. `VariantFactory`
Factory class for creating Variant objects from VCF records.

- **Responsibilities**:
  - Determine variant types (SNP, insertion, deletion, etc.)
  - Create appropriate Variant objects from VCF record data
  - Handle complex variant type detection logic

- **Key methods**:
  - `is_snp()`, `is_indel()`, `is_structural_variant()`: Static variant type detection methods
  - `get_variant_type()`: Determine the VariantType enum for a variant
  - `create_variants_from_record()`: Create Variant objects from a VcfRecord

### 7. `VcfVariant`
Higher-level representation of a variant.

- **Responsibilities**:
  - Provide a more user-friendly interface to variant data
  - Convert between VCF coordinates (1-based) and internal coordinates (0-based)
  - Represent structural variants, SNPs, indels, etc.

- **Key methods**:
  - `get_alleles()`: Get reference and alternate alleles
  - `get_genotypes()`: Get genotypes for all samples

### 8. `VcfIndex`
Provides indexing for random access to VCF files.

- **Responsibilities**:
  - Build or load index for VCF file
  - Support querying by genomic region

- **Key methods**:
  - `build_index()`: Create index for VCF file
  - `query(chrom, start, end)`: Get file positions for records in region

## Lazy Parsing Implementation

The key to efficient VCF parsing is to avoid parsing fields that aren't needed. Here's how lazy parsing will be implemented:

1. **Record Level Laziness**:
   - When iterating through a VCF file, only parse the fixed fields (CHROM, POS, etc.)
   - Store the raw string for INFO and FORMAT fields
   - Only parse these fields when explicitly requested

2. **INFO Field Laziness**:
   - Store the raw INFO string in the VcfRecord
   - When `get_info(field_id)` is called, parse only that specific field
   - Cache parsed results to avoid re-parsing

3. **FORMAT Field Laziness**:
   - Store the raw FORMAT string and sample data in the VcfRecord
   - When `get_format(field_id, sample_idx)` is called, parse only that specific field for that sample
   - Cache parsed results to avoid re-parsing

4. **Genotype Laziness**:
   - Only parse genotype data (GT field) when explicitly requested
   - Cache parsed genotypes

## Example Usage

```python
# Open a VCF file
vcf_reader = VcfReader("example.vcf.gz")

# Iterate through records
for record in vcf_reader:
    # Access fixed fields (parsed immediately)
    chrom = record.get_chrom()
    pos = record.get_pos()  # 0-based position
    
    # Lazy parsing of INFO fields - only parses AF when requested
    if record.has_info("AF"):
        allele_freq = record.get_info("AF")
    
    # Lazy parsing of genotypes - only parses when requested
    genotypes = record.get_genotypes()
    
    # Lazy parsing of FORMAT fields - only parses DP for sample 0 when requested
    if record.has_format("DP"):
        depth = record.get_format("DP", 0)
```

## Performance Considerations

1. **File Reading**:
   - Support both plain text and gzipped VCF files
   - Use buffered reading for better performance

2. **Memory Management**:
   - Avoid loading the entire file into memory
   - Use generators for iteration
   - Implement record caching for frequently accessed records

3. **Parsing Efficiency**:
   - Use specialized parsers for different field types
   - Cache parsed results to avoid re-parsing
   - Use string splitting and regular expressions efficiently

4. **Indexing**:
   - Support tabix indexing for compressed VCF files
   - Implement simple indexing for uncompressed files

This design provides a flexible and efficient framework for working with VCF files, with lazy parsing to minimize unnecessary computation and memory usage.