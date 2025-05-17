# FASTQ format

**FASTQ format** is a text-based [format](https://en.wikipedia.org/wiki/File_format "File format") for storing both a biological sequence (usually [nucleotide sequence](https://en.wikipedia.org/wiki/Nucleotide_sequence "Nucleotide sequence")) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single [ASCII](https://en.wikipedia.org/wiki/ASCII "ASCII") character for brevity.

It was originally developed at the [Wellcome Trust Sanger Institute](https://en.wikipedia.org/wiki/Wellcome_Trust_Sanger_Institute "Wellcome Trust Sanger Institute") to bundle a [FASTA formatted](https://en.wikipedia.org/wiki/FASTA_format "FASTA format") sequence and its quality data, but has become the _[de facto](https://en.wikipedia.org/wiki/De_facto "De facto")_ standard for storing the output of high-throughput sequencing instruments such as the [Illumina](https://en.wikipedia.org/wiki/Illumina_(company) "Illumina (company)") Genome Analyzer.[\[1\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Cock2009-1)

## Format

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=1 "Edit section: Format")\]

A FASTQ file has four line-separated fields per sequence:

- Field 1 begins with a '@' character and is followed by a sequence identifier and an _optional_ description (like a [FASTA](https://en.wikipedia.org/wiki/FASTA_format "FASTA format") title line).
- Field 2 is the raw sequence letters.
- Field 3 begins with a '+' character and is _optionally_ followed by the same sequence identifier (and any description) again.
- Field 4 encodes the quality values for the sequence in Field 2, and must contain the same number of symbols as letters in the sequence.

A FASTQ file containing a single sequence might look like this:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

```

The byte representing quality runs from 0x21 (lowest quality; '!' in ASCII) to 0x7e (highest quality; '~' in ASCII).
Here are the quality value characters in left-to-right increasing order of quality ( [ASCII](https://en.wikipedia.org/wiki/ASCII "ASCII")):

```
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

```

The original Sanger FASTQ files split long sequences and quality strings over multiple lines, as is typically done for [FASTA](https://en.wikipedia.org/wiki/FASTA_format "FASTA format") files. Accounting for this makes parsing more complicated due to the choice of "@" and "+" as markers (as these characters can also occur in the quality string). Multi-line FASTQ files (and consequently multi-line FASTQ parsers) are less common now that the majority of sequencing carried out is short-read [Illumina sequencing](https://en.wikipedia.org/wiki/Illumina_dye_sequencing "Illumina dye sequencing"), with typical sequence lengths of around 100bp.

### Illumina sequence identifiers

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=2 "Edit section: Illumina sequence identifiers")\]

Sequences from the [Illumina](https://en.wikipedia.org/wiki/Solexa "Solexa") software use a systematic identifier:

```
@HWUSI-EAS100R:6:73:941:1973#0/1

```

| HWUSI-EAS100R | the unique instrument name |
| 6 | flowcell lane |
| 73 | tile number within the flowcell lane |
| 941 | 'x'-coordinate of the cluster within the tile |
| 1973 | 'y'-coordinate of the cluster within the tile |
| #0 | index number for a multiplexed sample (0 for no indexing) |
| /1 | the member of a pair, /1 or /2 _(paired-end or mate-pair reads only)_ |

Versions of the Illumina pipeline since 1.4 appear to use **#NNNNNN** instead of **#0** for the multiplex ID, where **NNNNNN** is the sequence of the multiplex tag.

With Casava 1.8 the format of the '@' line has changed:

```
@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

```

| EAS139 | the unique instrument name |
| 136 | the run id |
| FC706VJ | the flowcell id |
| 2 | flowcell lane |
| 2104 | tile number within the flowcell lane |
| 15343 | 'x'-coordinate of the cluster within the tile |
| 197393 | 'y'-coordinate of the cluster within the tile |
| 1 | the member of a pair, 1 or 2 _(paired-end or mate-pair reads only)_ |
| Y | Y if the read is filtered (did not pass), N otherwise |
| 18 | 0 when none of the control bits are on, otherwise it is an even number |
| ATCACG | index sequence |

Note that more recent versions of Illumina software output a sample number (defined by the order of the samples in the sample sheet) in place of an index sequence when an index sequence is not explicitly specified for a sample in the sample sheet. For example, the following header might appear in a FASTQ file belonging to the first sample of a batch of samples:

```
@EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:1

```

### NCBI Sequence Read Archive

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=3 "Edit section: NCBI Sequence Read Archive")\]

FASTQ files from the [INSDC](https://en.wikipedia.org/wiki/International_Nucleotide_Sequence_Database_Collaboration "International Nucleotide Sequence Database Collaboration") [Sequence Read Archive](https://en.wikipedia.org/wiki/Sequence_Read_Archive "Sequence Read Archive") often include a description, e.g.

```
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC

```

In this example there is an NCBI-assigned identifier, and the description holds the original identifier from [Solexa/Illumina](https://en.wikipedia.org/wiki/Solexa "Solexa") (as described above) plus the read length. Sequencing was performed in paired-end mode (~500bp insert size), see [SRR001666](https://www.ncbi.nlm.nih.gov/sra/SRR001666). The default output format of fastq-dump produces entire spots, containing any technical reads and typically single or paired-end biological reads.

```
$ fastq-dump.2.9.0 -Z -X 2 SRR001666
Read 2 spots for SRR001666
Written 2 spots for SRR001666
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIIIIIII>IIIIII/
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=72
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=72
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I

```

Modern usage of FASTQ almost always involves splitting the spot into its biological reads, as described in submitter-provided metadata:

```
$ fastq-dump -X 2 SRR001666 --split-3
Read 2 spots for SRR001666
Written 2 spots for SRR001666
$ head SRR001666_1.fastq  SRR001666_2.fastq
==> SRR001666_1.fastq <==
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=36
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGA
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBI

==> SRR001666_2.fastq <==
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
AAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIDIIIIIII>IIIIII/
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=36
AGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT
+SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=36
IIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I

```

When present in the archive, fastq-dump can attempt to restore read names to original format. NCBI does not store original read names by default:

```
$ fastq-dump -X 2 SRR001666 --split-3 --origfmt
Read 2 spots for SRR001666
Written 2 spots for SRR001666
$ head SRR001666_1.fastq  SRR001666_2.fastq
==> SRR001666_1.fastq <==
@071112_SLXA-EAS1_s_7:5:1:817:345
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+071112_SLXA-EAS1_s_7:5:1:817:345
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
@071112_SLXA-EAS1_s_7:5:1:801:338
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGA
+071112_SLXA-EAS1_s_7:5:1:801:338
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBI

==> SRR001666_2.fastq <==
@071112_SLXA-EAS1_s_7:5:1:817:345
AAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
+071112_SLXA-EAS1_s_7:5:1:817:345
IIIIIIIIIIIIIIIIIIIIDIIIIIII>IIIIII/
@071112_SLXA-EAS1_s_7:5:1:801:338
AGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT
+071112_SLXA-EAS1_s_7:5:1:801:338
IIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I

```

In the example above, the original read names were used rather than the accessioned read name. NCBI accessions runs and the reads they contain. Original read names, assigned by sequencers, are able to function as locally unique identifiers of a read, and convey exactly as much information as a serial number. The ids above were algorithmically assigned based upon run information and geometric coordinates. Early SRA loaders parsed these ids and stored their decomposed components internally. NCBI stopped recording read names because they are frequently modified from the vendors' original format in order to associate some additional information meaningful to a particular processing pipeline, and this caused name format violations that resulted in a high number of rejected submissions. Without a clear schema for read names, their function remains that of a unique read id, conveying the same amount of information as a read serial number. See various [SRA Toolkit issues](https://github.com/ncbi/sra-tools/issues?q=is%3Aissue+original+read+names) for details and discussions.

Also note that [fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) converts this FASTQ data from the original Solexa/Illumina encoding to the Sanger standard (see encodings below). This is because [the SRA serves as a repository for NGS information, rather than format](https://github.com/ncbi/sra-tools/issues/130#issuecomment-409995254). The various \*-dump tools are capable of producing data in several formats from the same source. The requirements for doing so have been dictated by users over several years, with the majority of early demand coming from the [1000 Genomes Project](https://en.wikipedia.org/wiki/1000_Genomes_Project "1000 Genomes Project").

## Variations

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=4 "Edit section: Variations")\]

### Quality

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=5 "Edit section: Quality")\]

A quality value _Q_ is an integer mapping of _p_ (i.e., the probability that the corresponding base call is incorrect). Two different equations have been in use. The first is the standard Sanger variant to assess reliability of a base call, otherwise known as [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score "Phred quality score"):

Qsanger=−10log10⁡p{\\displaystyle Q\_{\\text{sanger}}=-10\\,\\log \_{10}p}![{\displaystyle Q_{\text{sanger}}=-10\,\log _{10}p}](https://wikimedia.org/api/rest_v1/media/math/render/svg/d35e9473c8fb52974f7f37c6bd66852e2e276da3)

The Solexa pipeline (i.e., the software delivered with the Illumina Genome Analyzer) earlier used a different mapping, encoding the [odds](https://en.wikipedia.org/wiki/Odds "Odds") _p_/(1- _p_) instead of the probability _p_:

Qsolexa-prior to v.1.3=−10log10⁡p1−p{\\displaystyle Q\_{\\text{solexa-prior to v.1.3}}=-10\\,\\log \_{10}{\\frac {p}{1-p}}}![{\displaystyle Q_{\text{solexa-prior to v.1.3}}=-10\,\log _{10}{\frac {p}{1-p}}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/ca73445cc28b13f2bde5997fe4499011de514bcc)

Although both mappings are asymptotically identical at higher quality values, they differ at lower quality levels (i.e., approximately _p_ \> 0.05, or equivalently, _Q_ < 13).

[![Relationship between Q and p](https://upload.wikimedia.org/wikipedia/commons/thumb/6/6b/Probability_metrics.svg/1200px-Probability_metrics.svg.png)](https://en.wikipedia.org/wiki/File:Probability_metrics.svg) Relationship between _Q_ and _p_ using the Sanger (red) and Solexa (black) equations (described above). The vertical dotted line indicates _p_ = 0.05, or equivalently, _Q_ ≈ 13.

At times there has been disagreement about which mapping Illumina actually uses. The user guide (Appendix B, page 122) for version 1.4 of the Illumina pipeline states that: "The scores are defined as ⁠Q=10⋅log10⁡p1−p{\\displaystyle Q=10\\cdot \\log \_{10}{\\tfrac {p}{1-p}}}![{\displaystyle Q=10\cdot \log _{10}{\tfrac {p}{1-p}}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/20fd051f508556b6e1d72adfd03903c8aa7a566d)⁠ \[ _[sic](https://en.wikipedia.org/wiki/Sic "Sic")_\], where p is the probability of a base call corresponding to the base in question".[\[2\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Illumina_User_Guide_1.4-2) In retrospect, this entry in the manual appears to have been an error. The user guide (What's New, page 5) for version 1.5 of the Illumina pipeline lists this description instead: "Important Changes in Pipeline v1.3 \[ _[sic](https://en.wikipedia.org/wiki/Sic "Sic")_\]. The quality scoring scheme has changed to the Phred \[i.e., Sanger\] scoring scheme, encoded as an ASCII character by adding 64 to the Phred value. A Phred score of a base is: Qphred=−10log10⁡e{\\displaystyle Q\_{\\text{phred}}=-10\\log \_{10}e}![{\displaystyle Q_{\text{phred}}=-10\log _{10}e}](https://wikimedia.org/api/rest_v1/media/math/render/svg/52711c8f2bbf23f584d917e086aaef72fde6d01c), where _e_ is the estimated probability of a base being wrong.[\[3\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Illumina_User_Guide_1.5-3)

### Encoding

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=6 "Edit section: Encoding")\]

- Sanger format can encode a [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score "Phred quality score") from 0 to 93 using ASCII 33 to 126 (although in raw read data the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps). Also used in SAM format.[\[4\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Sequence/Alignment_Map_format-4) Coming to the end of February 2011, Illumina's newest version (1.8) of their pipeline CASAVA will directly produce fastq in Sanger format, according to the announcement on seqanswers.com forum.[\[5\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Upcoming_changes_in_CASAVA_topic-5)
- Element Biosciences AVITI reads are encoded following the Sanger convention: Phred quality scores from 0 to 93 are encoded using ASCII 33 to 126. Raw reads typically exhibit base quality scores in the range of \[0, 55\]. [\[6\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Elembio_AVITI_FASTQ_format_specification-6)
- PacBio HiFi reads, which are typically stored in SAM/BAM format, use the Sanger convention: Phred quality scores from 0 to 93 are encoded using ASCII 33 to 126. Raw PacBio subreads use the same convention but typically assign a placeholder base quality (Q0) to all bases in the read.[\[7\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-PacBio_BAM_format_specification-7)
- Oxford Nanopore Duplex reads, called using the dorado basecaller are typically stored in SAM/BAM format. After changing to a 16-bit internal quality representation, the reported base quality limit is q50 (S).[\[8\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Duplex-dorado-8)
- Solexa/Illumina 1.0 format can encode a Solexa/Illumina quality score from -5 to 62 using [ASCII](https://en.wikipedia.org/wiki/ASCII "ASCII") 59 to 126 (although in raw read data Solexa scores from -5 to 40 only are expected)
- Starting with Illumina 1.3 and before Illumina 1.8, the format encoded a [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score "Phred quality score") from 0 to 62 using [ASCII](https://en.wikipedia.org/wiki/ASCII "ASCII") 64 to 126 (although in raw read data Phred scores from 0 to 40 only are expected).
- Starting in Illumina 1.5 and before Illumina 1.8, the Phred scores 0 to 2 have a slightly different meaning. The values 0 and 1 are no longer used and the value 2, encoded by ASCII 66 "B", is used also at the end of reads as a _Read Segment Quality Control Indicator_.[\[9\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-9) The Illumina manual[\[10\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-10) (page 30) states the following: _If a read ends with a segment of mostly low quality (Q15 or below), then all of the quality values in the segment are replaced with a value of 2 (encoded as the letter B in Illumina's text-based encoding of quality scores)... This Q2 indicator does not predict a specific error rate, but rather indicates that a specific final portion of the read should not be used in further analyses._ Also, the quality score encoded as "B" letter may occur internally within reads at least as late as pipeline version 1.6, as shown in the following example:

```
@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCTTGAGATTTGTTGGGGGAGACATTTTTGTGATTGCCTTGAT
+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]dddd`ddd^dddadd^BBBBBBBBBBBBBBBBBBBBBBBB

```

An alternative interpretation of this ASCII encoding has been proposed.[\[11\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-11) Also, in Illumina runs using PhiX controls, the character 'B' was observed to represent an "unknown quality score". The error rate of 'B' reads was roughly 3 phred scores lower the mean observed score of a given run.

- Starting in Illumina 1.8, the quality scores have basically returned to the use of the Sanger format (Phred+33).

For raw reads, the range of scores will depend on the technology and the base caller used, but will typically be up to 41 for recent Illumina chemistry. Since the maximum observed quality score was previously only 40, various scripts and tools break when they encounter data with quality values larger than 40. For processed reads, scores may be even higher. For example, quality values of 45 are observed in reads from Illumina's Long Read Sequencing Service (previously Moleculo).

```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...........................................
  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |              |               |                     |
 33                        59   64       73             88             104                   126
  0........................26...31.......40
                           -5....0........9.............................40
                                 0........9.............................40
                                    3.....9..............................41
  0.2......................26...31........41
  0..................20........30........40........50
  0..................20........30........40........50...55
  0..................20........30........40........50..........................................93

```

```
 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 41)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 N - Nanopore      Phred+33,  Duplex reads typically (0, 50)
 E - ElemBio AVITI Phred+33,  raw reads typically (0, 55)
 P - PacBio        Phred+33,  HiFi reads typically (0, 93)

```

### Color space

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=7 "Edit section: Color space")\]

For SOLiD data, the format is modified to a color space FASTQ sequence (CSFASTQ), where bases in the sequence are combined with the numbers 0, 1, 2, and 3, indicating how bases are modified relative to the previous base in the sequence (0: no change; 1: transition; 2: non-complementary transversion; 3: complementary transversion).[\[1\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Cock2009-1) This format matched the different sequencing chemistry used by SOLiD sequencers. Initial representations only used nucleotide bases at the start of the sequence, but later versions included bases embedded at periodic intervals to improve basecalling and mapping accuracy.

The quality values for CSFASTQ are identical to those of the Sanger format. Alignment tools differ in their preferred version of the quality values: some include a quality score (set to 0, i.e. '!') for the leading nucleotide, others do not. The sequence read archive includes this quality score.

### FAST5 and HDF5 evolutions

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=8 "Edit section: FAST5 and HDF5 evolutions")\]

The FAST4 format was invented as a derivative of the FASTQ format where each of the 4 bases (A,C,G,T) had separate probabilities stored. It was part of the [Swift](https://sourceforge.net/p/swiftng "sourceforge:p/swiftng") basecaller, an open source package for primary data analysis on next-gen sequence data "from images to basecalls".

The FAST5 format was invented as an extension of the FAST4 format. The FAST5 files are [Hierarchical Data Format 5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format "Hierarchical Data Format") (HDF5) files with a specific schema defined by [Oxford Nanopore Technologies](https://en.wikipedia.org/wiki/Oxford_Nanopore_Technologies "Oxford Nanopore Technologies") (ONT).[\[12\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-12)

### Simulation

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=9 "Edit section: Simulation")\]

FASTQ read simulation has been approached by several tools.[\[13\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-13)[\[14\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-14)
A comparison of those tools can be seen here.[\[15\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-15)

### Compression

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=10 "Edit section: Compression")\]

#### General compressors

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=11 "Edit section: General compressors")\]

General-purpose tools such as Gzip and bzip2 regard FASTQ as a plain text file and result in suboptimal compression ratios. NCBI's [Sequence Read Archive](https://en.wikipedia.org/wiki/Sequence_Read_Archive "Sequence Read Archive") encodes metadata using the LZ-77 scheme.
General FASTQ compressors typically compress distinct fields (read names, sequences, comments, and quality scores) in a FASTQ file separately; these include DSRC and DSRC2, FQC, LFQC, Fqzcomp, and Slimfastq.

#### Reads

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=12 "Edit section: Reads")\]

Having a reference genome around is convenient because then instead of storing the nucleotide sequences themselves, one can just align the reads to the reference genome and store the positions (pointers) and mismatches; the pointers can then be sorted according to their order in the reference sequence and encoded, e.g., with run-length encoding. When the [coverage](https://en.wikipedia.org/wiki/Coverage_(genetics) "Coverage (genetics)") or the repeat content of the sequenced genome is high, this leads to a high compression ratio.
Unlike the [SAM](https://en.wikipedia.org/wiki/SAM_(file_format) "SAM (file format)")/BAM formats, FASTQ files do not specify a reference genome. **Alignment-based FASTQ compressors** supports the use of either user-provided or _de novo_ assembled reference: LW-FQZip uses a provided reference genome and Quip, Leon, k-Path and KIC perform _**de novo**_ assembly using a [de Bruijn graph](https://en.wikipedia.org/wiki/De_Bruijn_graph "De Bruijn graph")-based approach.

Explicit read mapping and _de novo_ assembly are typically slow. **Reordering-based FASTQ compressors** first cluster reads that share long substrings and then independently compress reads in each cluster after reordering them or assembling them into longer [contigs](https://en.wikipedia.org/wiki/Contig "Contig"), achieving perhaps the best trade-off between the running time and compression rate. SCALCE is the first such tool, followed by Orcom and Mince. BEETL uses a generalized [Burrows–Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform "Burrows–Wheeler transform") for reordering reads, and HARC achieves better performance with hash-based reordering. AssemblTrie instead assembles reads into reference trees with as few total number of symbols as possible in the reference.[\[16\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-pmid29422526-16)[\[17\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Zhu_Numanagi%C4%87_Sahinalp_2018_pp._779-783-17)

Benchmarks for these tools are available.[\[18\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Numanagi%C4%87_Bonfield_Hach_Voges_pp._1005%E2%80%931008-18)

#### Quality values

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=13 "Edit section: Quality values")\]

Quality values account for about half of the required disk space in the FASTQ format (before compression), and therefore the compression of the quality values can significantly reduce storage requirements and speed up analysis and transmission of sequencing data. Both lossless and lossy compression are recently being considered in the literature. For example, the algorithm QualComp[\[19\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-19) performs lossy compression with a rate (number of bits per quality value) specified by the user. Based on rate-distortion theory results, it allocates the number of bits so as to minimize the MSE (mean squared error) between the original (uncompressed) and the reconstructed (after compression) quality values. Other algorithms for compression of quality values include SCALCE[\[20\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-20) and Fastqz.[\[21\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-21) Both are lossless compression algorithms that provide an optional controlled lossy transformation approach. For example, SCALCE reduces the alphabet size based on the observation that “neighboring” quality values are similar in general. For a benchmark, see.[\[22\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-Morteza-22)

As of the HiSeq 2500 Illumina gives the option to output qualities that have been coarse grained into quality bins. The binned scores are computed directly from the empirical quality score table, which is itself tied to the hardware, software and chemistry that were used during the sequencing experiment.[\[23\]](https://en.wikipedia.org/wiki/FASTQ_format#cite_note-23)

## File extension

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=14 "Edit section: File extension")\]

There is no standard [file extension](https://en.wikipedia.org/wiki/File_extension "File extension") for a FASTQ file, but .fq and .fastq are commonly used.

## Format converters

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=15 "Edit section: Format converters")\]

- [Biopython](https://en.wikipedia.org/wiki/Biopython "Biopython") version 1.51 onwards (interconverts Sanger, Solexa and Illumina 1.3+)
- [EMBOSS](https://en.wikipedia.org/wiki/EMBOSS "EMBOSS") version 6.1.0 patch 1 onwards (interconverts Sanger, Solexa and Illumina 1.3+)
- [BioPerl](https://en.wikipedia.org/wiki/BioPerl "BioPerl") version 1.6.1 onwards (interconverts Sanger, Solexa and Illumina 1.3+)
- [BioRuby](https://en.wikipedia.org/wiki/BioRuby "BioRuby") version 1.4.0 onwards (interconverts Sanger, Solexa and Illumina 1.3+)
- [BioJava](https://en.wikipedia.org/wiki/BioJava "BioJava") version 1.7.1 onwards (interconverts Sanger, Solexa and Illumina 1.3+)

## See also

\[ [edit](https://en.wikipedia.org/w/index.php?title=FASTQ_format&action=edit&section=16 "Edit section: See also")\]

- The [FASTA](https://en.wikipedia.org/wiki/FASTA_format "FASTA format") format, used to represent genome sequences.
- The [SAM](https://en.wikipedia.org/wiki/SAM_(file_format) "SAM (file format)") and [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format) "CRAM (file format)") formats, used to represent genome sequencer reads that have been aligned to genome sequences.
- The [GVF](https://en.wikipedia.org/w/index.php?title=GVF_format&action=edit&redlink=1 "GVF format (page does not exist)") format (Genome Variation Format), an extension based on the [GFF3](https://en.wikipedia.org/wiki/GFF3 "GFF3") format.
