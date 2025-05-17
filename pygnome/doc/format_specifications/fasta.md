# FASTA File format for DNA or protein sequences

| [Filename extensions](https://en.wikipedia.org/wiki/Filename_extension "Filename extension") | .fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn |
| [Internet media type](https://en.wikipedia.org/wiki/Media_type "Media type") | `text/x-fasta` |
| [Uniform Type Identifier (UTI)](https://en.wikipedia.org/wiki/Uniform_Type_Identifier "Uniform Type Identifier") | no |
| Developed by | [David J. Lipman](https://en.wikipedia.org/wiki/David_J._Lipman "David J. Lipman")<br>[William R. Pearson](https://en.wikipedia.org/wiki/William_Pearson_(scientist) "William Pearson (scientist)")[\[1\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-rapid-1)[\[2\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-improved-2) |
| Initial release | 1985 |
| Type of format | [Bioinformatics](https://en.wikipedia.org/wiki/Bioinformatics "Bioinformatics") |
| Extended from | [ASCII](https://en.wikipedia.org/wiki/ASCII "ASCII") for [FASTA](https://en.wikipedia.org/wiki/FASTA "FASTA") |
| Extended to | [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format "FASTQ format")[\[3\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-fastq-3) |
| Website | [www.ncbi.nlm.nih.gov/BLAST/fasta.shtml](https://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml) |

FASTA format

In [bioinformatics](https://en.wikipedia.org/wiki/Bioinformatics "Bioinformatics") and [biochemistry](https://en.wikipedia.org/wiki/Biochemistry "Biochemistry"), the **FASTA format** is a text-based [format](https://en.wikipedia.org/wiki/File_format "File format") for representing either [nucleotide sequences](https://en.wikipedia.org/wiki/Nucleotide_sequence "Nucleotide sequence") or amino acid (protein) sequences, in which nucleotides or [amino acids](https://en.wikipedia.org/wiki/Amino_acid "Amino acid") are represented using single-letter codes.

The format allows for sequence names and comments to precede the sequences. It originated from the [FASTA](https://en.wikipedia.org/wiki/FASTA "FASTA") software package and has since become a near-universal standard in [bioinformatics](https://en.wikipedia.org/wiki/Bioinformatics "Bioinformatics").[\[4\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-4)

The simplicity of FASTA format makes it easy to manipulate and parse sequences using text-processing tools and [scripting languages](https://en.wikipedia.org/wiki/Scripting_language "Scripting language").

## Overview

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=1 "Edit section: Overview")\]

A sequence begins with a greater-than character (">") followed by a description of the sequence (all in a single line). The lines immediately following the description line are the sequence representation, with one letter per amino acid or nucleic acid, and are typically no more than 80 characters in length.

For example:

```
>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*

```

### Original format

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=2 "Edit section: Original format")\]

The original FASTA/ [Pearson](https://en.wikipedia.org/wiki/William_Pearson_(scientist) "William Pearson (scientist)") format is described in the documentation for the [FASTA](https://en.wikipedia.org/wiki/FASTA "FASTA") suite of programs. It can be downloaded with any free distribution of FASTA (see fasta20.doc, fastaVN.doc, or fastaVN.me—where VN is the Version Number).

In the original format, a sequence was represented as a series of lines, each of which was no longer than 120 characters and usually did not exceed 80 characters. This probably was to allow for the preallocation of fixed line sizes in software: at the time most users relied on [Digital Equipment Corporation](https://en.wikipedia.org/wiki/Digital_Equipment_Corporation "Digital Equipment Corporation") (DEC) [VT220](https://en.wikipedia.org/wiki/VT220 "VT220") (or compatible) terminals which could display 80 or 132 characters per line.[\[5\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-5)[\[6\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-6) Most people preferred the bigger font in 80-character modes and so it became the recommended fashion to use 80 characters or less (often 70) in FASTA lines. Also, the width of a standard printed page is 70 to 80 characters (depending on the font). Hence, 80 characters became the norm.[\[7\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-7)

The first line in a FASTA file started either with a ">" (greater-than) symbol or, less frequently, a ";"[\[8\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-8) (semicolon) was taken as a comment. Subsequent lines starting with a semicolon would be ignored by software. Since the only comment used was the first, it quickly became used to hold a summary description of the sequence, often starting with a unique library accession number, and with time it has become commonplace to always use ">" for the first line and to not use ";" comments (which would otherwise be ignored).

Following the initial line (used for a unique description of the sequence) was the actual sequence itself in the standard one-letter character string. Anything other than a valid character would be ignored (including spaces, tabulators, asterisks, etc...). It was also common to end the sequence with an "\*" (asterisk) character (in analogy with use in PIR formatted sequences) and, for the same reason, to leave a blank line between the description and the sequence. Below are a few sample sequences:

```
;LCBO - Prolactin precursor - Bovine
; a sample sequence in FASTA format
MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSS
EMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL
VTEVRGMKGAPDAILSRAIEIEEENKRLLEGMEMIFGQVIPGAKETEPYPVWSGLPSLQTKDED
ARYSAFYNLLHCLRRDSSKIDTYLKLLNCRIIYNNNC*

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*

>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

```

A multiple-sequence FASTA format, or multi-FASTA format, would be obtained by concatenating several single-sequence FASTA files in one file. This does not imply a contradiction with the format as only the first line in a FASTA file may start with a ";" or ">", forcing all subsequent sequences to start with a ">" in order to be taken as separate sequences (and further forcing the exclusive reservation of ">" for the sequence definition line). Thus, the examples above would be a multi-FASTA file if taken together.

Modern bioinformatics programs that rely on the FASTA format expect the sequence headers to be preceded by ">". The sequence is generally represented as "interleaved", or on multiple lines as in the above example, but may also be "sequential", or on a single line. Running different bioinformatics programs may require conversions between "sequential" and "interleaved" FASTA formats.

## Description line

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=3 "Edit section: Description line")\]

The description line (defline) or header/identifier line, which begins with ">", gives a name and/or a unique identifier for the sequence, and may also contain additional information. In a deprecated practice, the header line sometimes contained more than one header, separated by a ^A (Control-A) character. In the original [Pearson](https://en.wikipedia.org/wiki/William_Pearson_(scientist) "William Pearson (scientist)") FASTA format, one or more comments, distinguished by a semi-colon at the beginning of the line, may occur after the header. Some databases and bioinformatics applications do not recognize these comments and follow [the NCBI FASTA specification](https://www.ncbi.nlm.nih.gov/blast/fasta.shtml). An example of a multiple sequence FASTA file follows:

```
>SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH

```

### NCBI identifiers

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=4 "Edit section: NCBI identifiers")\]

The [NCBI](https://en.wikipedia.org/wiki/National_Center_for_Biotechnology_Information "National Center for Biotechnology Information") defined a standard for the unique identifier used for the sequence (SeqID) in the header line. This allows a sequence that was obtained from a database to be labelled with a reference to its database record. The database identifier format is understood by the NCBI tools like `makeblastdb` and `table2asn`. The following list describes the NCBI FASTA defined format for sequence identifiers.[\[9\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-9)

| Type | Format(s) | Example(s) |
| --- | --- | --- |
| local (i.e. no database reference) | `lcl|integer`<br>`lcl|string` | `lcl|123`<br>`lcl|hmm271` |
| GenInfo backbone seqid | `bbs|integer` | `bbs|123` |
| GenInfo backbone moltype | `bbm|integer` | `bbm|123` |
| GenInfo import ID | `gim|integer` | `gim|123` |
| [GenBank](https://www.ncbi.nlm.nih.gov/Genbank/index.html) | `gb|accession|locus` | `gb|M73307|AGMA13GT` |
| [EMBL](http://www.embl-heidelberg.de/) | `emb|accession|locus` | `emb|CAM43271.1|` |
| [PIR](https://web.archive.org/web/20140312021627/http://pir.georgetown.edu/) | `pir|accession|name` | `pir||G36364` |
| [SWISS-PROT](http://www.ebi.ac.uk/swissprot) | `sp|accession|name` | `sp|P01013|OVAX_CHICK` |
| patent | `pat|country|patent|sequence-number` | `pat|US|RE33188|1` |
| pre-grant patent | `pgp|country|application-number|sequence-number` | `pgp|EP|0238993|7` |
| [RefSeq](https://www.ncbi.nlm.nih.gov/projects/RefSeq) | `ref|accession|name` | `ref|NM_010450.1|` |
| general database reference<br>(a reference to a database that's not in this list) | `gnl|database|integer`<br>`gnl|database|string` | `gnl|taxon|9606`<br>`gnl|PID|e1632` |
| GenInfo integrated database | `gi|integer` | `gi|21434723` |
| [DDBJ](http://www.ddbj.nig.ac.jp/) | `dbj|accession|locus` | `dbj|BAC85684.1|` |
| [PRF](http://www.prf.or.jp/) | `prf|accession|name` | `prf||0806162C` |
| [PDB](https://web.archive.org/web/20080828002005/http://www.rcsb.org./pdb) | `pdb|entry|chain` | `pdb|1I4L|D` |
| third-party [GenBank](https://www.ncbi.nlm.nih.gov/Genbank/index.html) | `tpg|accession|name` | `tpg|BK003456|` |
| third-party [EMBL](http://www.embl-heidelberg.de/) | `tpe|accession|name` | `tpe|BN000123|` |
| third-party [DDBJ](http://www.ddbj.nig.ac.jp/) | `tpd|accession|name` | `tpd|FAA00017|` |
| TrEMBL | `tr|accession|name` | `tr|Q90RT2|Q90RT2_9HIV1` |

The vertical bars ("\|") in the above list are not separators in the sense of the [Backus–Naur form](https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form "Backus–Naur form") but are part of the format. Multiple identifiers can be concatenated, also separated by vertical bars.

## Sequence representation

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=5 "Edit section: Sequence representation")\]

Following the header line, the actual sequence is represented. Sequences may be [protein sequences](https://en.wikipedia.org/wiki/Primary_structure "Primary structure") or [nucleic acid](https://en.wikipedia.org/wiki/Nucleic_acid "Nucleic acid") sequences, and they can contain gaps or alignment characters (see [sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment "Sequence alignment")). Sequences are expected to be represented in the standard IUB/IUPAC [amino acid](https://en.wikipedia.org/wiki/Amino_acid "Amino acid") and [nucleic acid](https://en.wikipedia.org/wiki/Nucleic_acid "Nucleic acid") codes, with these exceptions: lower-case letters are accepted and are mapped into upper-case; a single hyphen or dash can be used to represent a gap character; and in amino acid sequences, U and \* are acceptable letters (see below). Numerical digits are not allowed but are used in some databases to indicate the position in the sequence. The nucleic acid codes supported are:[\[10\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-10)[\[11\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-11)[\[12\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-12)

| Nucleic Acid Code | Meaning | Mnemonic |
| --- | --- | --- |
| A | A | [**A** denine](https://en.wikipedia.org/wiki/Adenine "Adenine") |
| C | C | [**C** ytosine](https://en.wikipedia.org/wiki/Cytosine "Cytosine") |
| G | G | [**G** uanine](https://en.wikipedia.org/wiki/Guanine "Guanine") |
| T | T | [**T** hymine](https://en.wikipedia.org/wiki/Thymine "Thymine") |
| U | U | [**U** racil](https://en.wikipedia.org/wiki/Uracil "Uracil") |
| (i) | i | [**i** nosine](https://en.wikipedia.org/wiki/Inosine "Inosine") (non-standard) |
| R | A or G (I) | [pu **R** ine](https://en.wikipedia.org/wiki/Purine "Purine") |
| Y | C, T or U | [p **Y** rimidines](https://en.wikipedia.org/wiki/Pyrimidine "Pyrimidine") |
| K | G, T or U | bases which are [**K** etones](https://en.wikipedia.org/wiki/Ketone "Ketone") |
| M | A or C | bases with [a **M** ino groups](https://en.wikipedia.org/wiki/Amino "Amino") |
| S | C or G | **S** trong interaction |
| W | A, T or U | **W** eak interaction |
| B | not A (i.e. C, G, T or U) | **B** comes after A |
| D | not C (i.e. A, G, T or U) | **D** comes after C |
| H | not G (i.e., A, C, T or U) | **H** comes after G |
| V | neither T nor U (i.e. A, C or G) | **V** comes after U |
| N | A C G T U | **N** ucleic acid |
| - | gap of indeterminate length |  |

The amino acid codes supported (22 amino acids and 3 special codes) are:

| Amino Acid Code | Meaning |
| --- | --- |
| A | [Alanine](https://en.wikipedia.org/wiki/Alanine "Alanine") |
| B | [Aspartic acid](https://en.wikipedia.org/wiki/Aspartic_acid "Aspartic acid") (D) or [Asparagine](https://en.wikipedia.org/wiki/Asparagine "Asparagine") (N) |
| C | [Cysteine](https://en.wikipedia.org/wiki/Cysteine "Cysteine") |
| D | [Aspartic acid](https://en.wikipedia.org/wiki/Aspartic_acid "Aspartic acid") |
| E | [Glutamic acid](https://en.wikipedia.org/wiki/Glutamic_acid "Glutamic acid") |
| F | [Phenylalanine](https://en.wikipedia.org/wiki/Phenylalanine "Phenylalanine") |
| G | [Glycine](https://en.wikipedia.org/wiki/Glycine "Glycine") |
| H | [Histidine](https://en.wikipedia.org/wiki/Histidine "Histidine") |
| I | [Isoleucine](https://en.wikipedia.org/wiki/Isoleucine "Isoleucine") |
| J | [Leucine](https://en.wikipedia.org/wiki/Leucine "Leucine") (L) or [Isoleucine](https://en.wikipedia.org/wiki/Isoleucine "Isoleucine") (I) |
| K | [Lysine](https://en.wikipedia.org/wiki/Lysine "Lysine") |
| L | [Leucine](https://en.wikipedia.org/wiki/Leucine "Leucine") |
| M | [Methionine](https://en.wikipedia.org/wiki/Methionine "Methionine")/ [Start codon](https://en.wikipedia.org/wiki/Start_codon "Start codon") |
| N | [Asparagine](https://en.wikipedia.org/wiki/Asparagine "Asparagine") |
| O | [Pyrrolysine](https://en.wikipedia.org/wiki/Pyrrolysine "Pyrrolysine") (rare) |
| P | [Proline](https://en.wikipedia.org/wiki/Proline "Proline") |
| Q | [Glutamine](https://en.wikipedia.org/wiki/Glutamine "Glutamine") |
| R | [Arginine](https://en.wikipedia.org/wiki/Arginine "Arginine") |
| S | [Serine](https://en.wikipedia.org/wiki/Serine "Serine") |
| T | [Threonine](https://en.wikipedia.org/wiki/Threonine "Threonine") |
| U | [Selenocysteine](https://en.wikipedia.org/wiki/Selenocysteine "Selenocysteine") (rare) |
| V | [Valine](https://en.wikipedia.org/wiki/Valine "Valine") |
| W | [Tryptophan](https://en.wikipedia.org/wiki/Tryptophan "Tryptophan") |
| Y | [Tyrosine](https://en.wikipedia.org/wiki/Tyrosine "Tyrosine") |
| Z | [Glutamic acid](https://en.wikipedia.org/wiki/Glutamic_acid "Glutamic acid") (E) or [Glutamine](https://en.wikipedia.org/wiki/Glutamine "Glutamine") (Q) |
| X | any |
| \* | translation stop |
| - | gap of indeterminate length |

## FASTA file

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=6 "Edit section: FASTA file")\]

### Filename extension

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=7 "Edit section: Filename extension")\]

There is no standard [filename extension](https://en.wikipedia.org/wiki/Filename_extension "Filename extension") for a text file containing FASTA formatted sequences. The table below shows each extension and its respective meaning.

| Extension | Meaning | Notes |
| --- | --- | --- |
| fasta, fas, fa[\[13\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-13) | generic FASTA | Any generic FASTA file |
| fna | FASTA nucleic acid | Used generically to specify nucleic acids |
| ffn | FASTA nucleotide of gene regions | Contains coding regions for a genome |
| faa | FASTA amino acid | Contains amino acid sequences |
| mpfa | FASTA amino acids | Contains multiple protein sequences |
| frn | FASTA [non-coding RNA](https://en.wikipedia.org/wiki/Non-coding_RNA "Non-coding RNA") | Contains non-coding RNA regions for a genome, e.g. tRNA, rRNA |

### Compression

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=8 "Edit section: Compression")\]

The compression of FASTA files requires a specific compressor to handle both channels of information: identifiers and sequence. For improved compression results, these are mainly divided into two streams where the compression is made assuming independence. For example, the algorithm MFCompress[\[14\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-MFCompress-14) performs lossless compression of these files using context modelling and arithmetic encoding. Genozip,[\[15\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-Genozip-15) a software package for compressing genomic files, uses an extensible context-based model. Benchmarks of FASTA file compression algorithms have been reported by Hosseini et al. in 2016,[\[16\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-Morteza-16) and Kryukov et al. in 2020.[\[17\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-SCB-17)

### Encryption

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=9 "Edit section: Encryption")\]

The encryption of FASTA files can be performed with various tools, including Cryfa and Genozip. Cryfa uses AES encryption and also enables data compression.[\[18\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-CRYFA1-18)[\[19\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-CRYFA2-19) Similarly, Genozip can encrypt FASTA files with AES-256 during compression.[\[15\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-Genozip-15)

## Extensions

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=10 "Edit section: Extensions")\]

[FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format "FASTQ format") is a form of FASTA format extended to indicate information related to sequencing. It is created by the [Sanger Centre](https://en.wikipedia.org/wiki/Sanger_Centre "Sanger Centre") in Cambridge.[\[3\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-fastq-3)

A2M/A3M are a family of FASTA-derived formats used for [sequence alignments](https://en.wikipedia.org/wiki/Sequence_alignment "Sequence alignment"). In A2M/A3M sequences, lowercase characters are taken to mean insertions, which are then indicated in the other sequences as the dot (".") character. The dots can be discarded for compactness without loss of information. As with typical FASTA files used in alignments, the gap ("-") is taken to mean exactly one position.[\[20\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-20) A3M is similar to A2M, with the added rule that gaps aligned to insertions can too be discarded.[\[21\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-21)

## Working with FASTA files

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=11 "Edit section: Working with FASTA files")\]

A plethora of user-friendly scripts are available from the community to perform FASTA file manipulations. Online toolboxes, such as FaBox[\[22\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-FaBox-22) or the FASTX-Toolkit within Galaxy servers, are also available.[\[23\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-Galaxyserver-23) These can be used to segregate sequence headers/identifiers, rename them, shorten them, or extract sequences of interest from large FASTA files based on a list of wanted identifiers (among other available functions). A tree-based approach to sorting multi-FASTA files (TREE2FASTA[\[24\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-tree2fasta-24)) also exists based on the coloring and/or annotation of sequences of interest in the FigTree viewer. Additionally, the [Bioconductor](https://en.wikipedia.org/wiki/Bioconductor "Bioconductor") _Biostrings_ package can be used to read and manipulate FASTA files in [R](https://en.wikipedia.org/wiki/R_(programming_language) "R (programming language)").[\[25\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-25)

Several online format converters exist to rapidly reformat multi-FASTA files to different formats (e.g. NEXUS, PHYLIP) for use with different phylogenetic programs, such as the converter available on phylogeny.fr.[\[26\]](https://en.wikipedia.org/wiki/FASTA_format#cite_note-phylodotfr-26)

## See also

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTA_format&action=edit&section=12 "Edit section: See also")\]

- The [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format "FASTQ format"), used to represent DNA sequencer reads along with quality scores.
- The [SAM](https://en.wikipedia.org/wiki/SAM_(file_format) "SAM (file format)") and [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format) "CRAM (file format)") formats, used to represent genome sequencer reads that have been aligned to genome sequences.
- The GVF format (Genome Variation Format), an extension based on the [GFF3](https://en.wikipedia.org/wiki/GFF3 "GFF3") format.
