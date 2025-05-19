# GTF Parsing and Transcript Structure Test Plan

This test plan outlines the test cases for verifying that transcript structures are correctly created from GTF files, with a focus on the "lazy" creation of various transcript features (CDS, UTR, splice sites, etc.).

## Test Cases

### 1. Basic Transcript Structure Creation
- [x] **Forward Strand**: Test that a basic transcript with exons is correctly parsed from GTF
- [x] **Reverse Strand**: Test that a basic transcript with exons is correctly parsed from GTF on the reverse strand

### 2. CDS Creation and Handling
- [x] **Forward Strand with CDS**: Test that CDS regions are correctly parsed and added to the transcript
- [x] **Reverse Strand with CDS**: Test that CDS regions are correctly parsed and added to the transcript on the reverse strand
- [x] **No CDS Provided**: Test that a transcript without CDS regions is handled correctly
- [x] **CDS with Phase Information**: Test that phase information is correctly parsed and used
- [x] **CDS without Phase Information**: Test that phase is correctly inferred when not provided

### 3. UTR Creation and Handling
- [x] **Forward Strand with UTRs**: Test that 5' and 3' UTRs are correctly parsed and added to the transcript
- [x] **Reverse Strand with UTRs**: Test that 5' and 3' UTRs are correctly parsed and added to the transcript on the reverse strand
- [x] **No UTRs Provided**: Test that a transcript without UTRs is handled correctly
- [x] **CDS without UTRs**: Test that a transcript with CDS but no explicit UTRs correctly infers UTRs

### 4. Splice Site Creation
- [x] **Forward Strand Splice Sites**: Test that splice sites are correctly inferred from exon boundaries
- [x] **Reverse Strand Splice Sites**: Test that splice sites are correctly inferred from exon boundaries on the reverse strand
- [x] **Single Exon Transcript**: Test that a transcript with a single exon has no splice sites

### 5. Intron Creation
- [x] **Forward Strand Introns**: Test that introns are correctly inferred from exon boundaries
- [x] **Reverse Strand Introns**: Test that introns are correctly inferred from exon boundaries on the reverse strand
- [x] **Single Exon Transcript**: Test that a transcript with a single exon has no introns

### 6. Protein Translation
- [x] **Forward Strand Translation**: Test that a transcript with CDS correctly translates to protein
- [x] **Reverse Strand Translation**: Test that a transcript with CDS on the reverse strand correctly translates to protein
- [x] **Non-coding Transcript**: Test that a non-coding transcript returns an empty protein sequence

### 7. Edge Cases
- [ ] **Overlapping Exons**: Test that overlapping exons are handled correctly
- [ ] **Overlapping CDS**: Test that overlapping CDS regions are handled correctly
- [x] **CDS Spanning Multiple Exons**: Test that CDS regions spanning multiple exons are handled correctly
- [ ] **UTRs Spanning Multiple Exons**: Test that UTRs spanning multiple exons are handled correctly
- [ ] **Incomplete CDS (not multiple of 3)**: Test that CDS with length not divisible by 3 is handled correctly