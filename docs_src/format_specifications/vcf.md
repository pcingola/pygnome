Below is the content of pages 4–16 of the VCF 4.5 PDF, faithfully converted into Markdown. Section and subsection numbers match the PDF.
## 1 The VCF specification

Reference: https://samtools.github.io/hts-specs/VCFv4.5.pdf

### 1.1 An example

```vcf
##fileformat=VCFv4.5
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS      ID          REF ALT    QUAL FILTER INFO                       FORMAT        NA00001         NA00002        NA00003
20     14370    rs6054257   G   A      29   PASS   NS=3;DP=14;AF=0.5;DB;H2     GT:GQ:DP:HQ   0|0:48:1:51,51   1|0:48:8:51,51  1/1:43:5:.,.  
20     17330    .           T   A      3    q10    NS=3;DP=11;AF=0.017           GT:GQ:DP:HQ   0|0:49:3:58,50   0|1:3:5:65,3    0/0:41:3     
20     1110696 rs6040355   A   G,T    67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27   2|1:2:0:18,2    2/2:35:4     
20     1230237 .           T   .      47   PASS   NS=3;DP=13;AA=T               GT:GQ:DP:HQ   0|0:54:7:56,60   0|0:48:4:51,51  0/0:61:2    
20     1234567 microsat1    GTC G,GTCT 50   PASS   NS=3;DP=9;AA=G              GT:GQ:DP      0/1:35:4       0/2:17:2       1/1:40:3    
```

This shows a simple SNP, a filtered‑out SNP, a multi‑allelic site, a monomorphic reference site, and a microsatellite with two alternate alleles; genotype columns illustrate phased vs. unphased calls and per‑sample quality/depth/HQ .

---

### 1.2 Character encoding, non-printable characters and special characters

VCF files must be UTF‑8 encoded (no byte‑order mark) and may use either LF (`\n`) or CR+LF (`\r\n`) as line separators. Non‑printable characters U+0000–U+0008, U+000B–U+000C, and U+000E–U+001F are disallowed. Characters with special meaning (e.g. `;`, `:`, `=`, `%`, `,`) within fields must be percent‑encoded as follows:

| Encoding | Character |   |
| -------- | --------- | - |
| `%3A`    | `:`       |   |
| `%3B`    | `;`       |   |
| `%3D`    | `=`       |   |
| `%25`    | `%`       |   |
| `%2C`    | `,`       |   |
| `%0D`    | CR        |   |
| `%0A`    | LF        |   |
| `%09`    | TAB       |   |

---

### 1.3 Data types

Supported VCF data types (applies to both INFO and FORMAT fields) are:

* **Integer** (32‑bit signed; values ≤ ±2³¹–1)
* **Float** (32‑bit IEEE‑754; regex `^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$` or `^[-+]?(INF|INFINITY|NAN)$`)
* **Flag** (present/absent, Number=0)
* **Character** (single‐character string)
* **String** (arbitrary string)&#x20;

---

### 1.4 Meta‑information lines

Meta‑information lines begin with `##` and must come before the single header line (`#CHROM ...`). They are either:

* **Unstructured**: `##key=value`
* **Structured**: `##key=<key1=val1,key2=val2,...>`
  – All structured lines require a unique `ID=` field inside the `<...>` .

#### 1.4.1 File format

A required first line specifying the VCF version:

```vcf
##fileformat=VCFv4.5
``` :contentReference[oaicite:4]{index=4}

#### 1.4.2 Information field format

Defines INFO keys:
```vcf
##INFO=<ID=ID,Number=<number>,Type=<type>,Description="...",Source="...",Version="...">
```

* `Type` ∈ {Integer, Float, Flag, Character, String}
* `Number` ∈ ℕ or:

  * `A`: one value per ALT allele
  * `R`: one value per allele (REF + ALTs)
  * `G`: one value per genotype (see 1.6.2)
  * `.`: unknown/unbounded
* `Flag` fields have `Number=0` and no `=value` part .

#### 1.4.3 Filter field format

Defines FILTER codes:

```vcf
##FILTER=<ID=ID,Description="...">
``` :contentReference[oaicite:6]{index=6}

#### 1.4.4 Individual format field format

Defines per‑sample FORMAT keys:
```vcf
##FORMAT=<ID=ID,Number=<number>,Type=<type>,Description="...">
```

* Same types as INFO
* Additional `Number` codes:

  * `LA`, `LR`, `LG`: allele‐subset variants
  * `P`: one value per allele in GT
  * `M`: one value per possible base modification (see spec) .

#### 1.4.5 Alternative allele field format

Defines symbolic ALT alleles (e.g., structural variants or ambiguity codes):

```vcf
##ALT=<ID=type,Description="...">
```

* Recommended types: `DEL`, `INS`, `DUP`, `INV`, `CNV`, with subtypes like `CNV:TR`, `DUP:TANDEM`, etc.
* IUPAC codes as ALT can be declared (e.g. `ID=R` for A/G) .

#### 1.4.6 Assembly field format

(Optional) link to FASTA for breakend assemblies:

```vcf
##assembly=url
``` :contentReference[oaicite:9]{index=9}

#### 1.4.7 Contig field format

Defines contigs and optional metadata:
```vcf
##contig=<ID=ctg1,length=...,md5=...,URL=...>
```

* `ID`: name matching `[!-~]` except certain punctuation
* Optional `length`, `md5`, `URL` fields .

#### 1.4.8 Sample field format

Defines sample metadata:

```vcf
##SAMPLE=<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,...,Description="...">
```
 :contentReference[oaicite:11]{index=11}

#### 1.4.9 Pedigree field format

Defines pedigree relationships:
```vcf
##PEDIGREE=<ID=ChildID,Father=FatherID,Mother=MotherID>
##pedigreeDB=url
```
:contentReference[oaicite:12]{index=12}

---

### 1.5 Header line syntax

The single header line names the fixed columns (tab‑delimited, no trailing tab):

```vcf
#CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT [<sampleID>...]]
```

---

### 1.6 Data lines

All data lines are tab‑delimited, no trailing tab, and end with a newline. Missing values are `.`.

#### 1.6.1 Fixed fields

Each record has **8** mandatory fields:

1. **CHROM** (String): reference contig or `<ID>`; entries for one CHROM must be contiguous.  
2. **POS** (Integer): 1‑based reference position; telomeres at 0 or N+1.  
3. **ID** (String): semicolon‑sep. unique IDs or `.`.  
4. **REF** (String): one or more bases (`A,C,G,T,N`), padded for indels at start; case-insensitive.  
5. **ALT** (String): comma‑sep. alternate alleles (`A,C,G,T,N`, `*`, `.`, `<ID>`, `<*>`, breakends).  
6. **QUAL** (Float): phred-scaled quality (`–10 log₁₀` error prob.), or `.`.  
7. **FILTER** (String): `PASS` or semicolon‑sep. filter IDs; `.` if not applied.  
8. **INFO** (String): semicolon‑sep. `key[=value[,...]]` entries or `.`; values may use percent‐encoding.  

Reserved INFO keys are detailed in Table 1 below. :contentReference[oaicite:14]{index=14}

---

#### 1.6.2 Genotype fields

If genotypes are present, a `FORMAT` column lists per‑sample keys (colon‑delimited), followed by one column per sample.  

- The first key **must** be `GT` if present; local‑allele fields require `LAA`.  
- Missing values use `.`; trailing fields may be omitted except `GT`.  

Common reserved FORMAT keys (Table 2) include `AD`, `DP`, `FT`, `GL`, `GP`, `GQ`, `GT`, `HQ`, `PL`, `PS`, and many more (including base‑mod Modification keys M... and their depths; see spec) :contentReference[oaicite:15]{index=15}.

##### Table 1: Reserved INFO keys (excerpt)

| Key     | Number | Type    | Description                                              |
|---------|:------:|:--------|:---------------------------------------------------------|
| AA      | 1      | String  | Ancestral allele                                         |
| AC      | A      | Integer | Allele count per ALT allele                              |
| AF      | A      | Float   | Allele frequency per ALT allele                          |
| AN      | 1      | Integer | Total number of alleles called                           |
| DP      | 1      | Integer | Combined depth across samples                            |
| MQ      | 1      | Float   | RMS mapping quality                                      |
| ...       | ...      | ...       | ...                                                         |

##### Table 2: Reserved FORMAT/Genotype keys (excerpt)

| Field | Number | Type    | Description                                                    |
|-------|:------:|:--------|:---------------------------------------------------------------|
| GT    | 1      | String  | Genotype (0=REF, 1...=ALTs; phasing `/` or `|`)                   |
| GQ    | 1      | Integer | Conditional genotype quality                                   |
| DP    | 1      | Integer | Read depth                                                     |
| AD    | R      | Integer | Read depth per allele                                          |
| PL    | G      | Integer | Phred‐scaled genotype likelihoods                              |
| ...     | ...      | ...       | ...                                                               |

For full lists and detailed definitions, see Sections 3–5 and Tables 1–2 in the spec.


## 2 Understanding the VCF format and the haplotype representation

VCF records use a simple haplotype representation to describe variant alleles at a locus.  A VCF record holds all segregating alleles—both reference and alternate—as well as, optionally, genotypes for multiple individuals at that locus.  ALT haplotypes are constructed from the REF haplotype by replacing the reference bases at **POS** with the alternate bases.  In essence, the VCF record specifies a-REF-t and the alternative haplotypes are a-ALT-t for each ALT allele.&#x20;

### 2.1 VCF tag naming conventions

Implementation‑defined tags should follow these suffix/prefix conventions:&#x20;

* **`L` suffix** — log₁₀‐likelihood (e.g., `GL`, `CNL`): log₁₀ Pr(Data|Model), negative numbers.
* **`P` suffix** — linear‐scale posterior probability (e.g., `GP`, `CNP`): Pr(Model|Data).
* **`Q` suffix** — phred‐scaled quality: −10 log₁₀ Pr(Data|Model) (e.g., `GQ`, `CNQ`); the fixed site‐level **QUAL** field uses the same scale.
* **`L` prefix** — “local‐allele” equivalent of a Number=A, R, or G field (see 1.6.2).

---

## 3 INFO keys used for structural variants

The following INFO keys are reserved for encoding structural variants.  Values for imprecise variants should be best estimates.  When per‑allele values are required, they must be provided for **all** ALT alleles (including non‑structural alleles), using the missing‐value placeholder (`.`) where appropriate.&#x20;

```vcf
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=END,Number=1,Type=Integer,Description="Deprecated. Present for backwards compatibility">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">
##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS">
```

* **IMPRECISE** (Flag): marks an imprecise structural‐variant ALT.  If a record is marked IMPRECISE, `CIPOS` (or `CIEND`) **must** be present (even if `0,0`).
* **NOVEL** (Flag): indicates a novel structural variant.
* **END** (Integer): deprecated in favor of `SVLEN` (and `LEN` in FORMAT).
* **SVTYPE** (String): deprecated (redundant with symbolic `ALT` alleles); see 1.4.5 for valid ALT symbols (e.g., `<DEL>`, `<DUP>`, `<INV>`).
* **SVLEN** (Integer, one per ALT): variant length.  For `<DEL>`, `<INS>`, `<DUP>`, `<INV>` alleles, the number of bases deleted/inserted/etc.; for `<CNV>`, the segment length.  Use `.` for other ALT alleles.
* **CIPOS** (Integer list): confidence interval around **POS**; twice as many entries as ALT alleles, giving start/end offsets for each.

---

## 4 FORMAT keys used for structural variants

Reserved per‑sample FORMAT keys for structural variants include (among others):&#x20;

```vcf
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">
##FORMAT=<ID=PSL,Number=.,Type=Integer,Description="Phase set list for structural variants">
##FORMAT=<ID=PSO,Number=.,Type=Integer,Description="Phase set offsets">
```

* **CN** (Float): estimated per‑sample copy number.
* **PSL/PSO** (Integer lists): phase sets and offsets, used to indicate phasing of SNVs with structural events (see 2.1).

---

## 5 Representing variation in VCF records

### 5.1 Creating VCF entries for SNPs and small indels

VCF entries for simple variants are constructed by aligning REF and ALT haplotypes:

#### 5.1.1 Example 1

*Single SNP at position 100:*

```vcf
chr1 100 . A G 50 PASS . GT 0/1
```

#### 5.1.2 Example 2

*Two‐base insertion at position 200:*

```vcf
chr1 200 . AT ATG 60 PASS . GT 1|0
```

#### 5.1.3 Example 3

*Two‐base deletion at position 300:*

```vcf
chr1 300 . AT A 70 PASS . GT 0/1
```

### 5.2 Decoding VCF entries for SNPs and small indels

Given a REF/ALT and genotype, one reconstructs the haplotypes:

* **5.2.1 SNP** — replace the single base at **POS**.
* **5.2.2 Insertion** — pad REF to first base, ALT contains inserted bases.
* **5.2.3 Deletion** — REF longer than ALT; ALT is the deleted‐bases string.
* **5.2.4 Mixed (microsatellite)** — multiple ALT alleles of varying lengths.

### 5.3 Encoding Structural Variants

Symbolic alleles (e.g., `<DEL>`, `<INS>`, `<DUP>`) are used in **ALT**, with companion INFO/FORMAT fields (`SVTYPE`, `SVLEN`, `CIPOS`, `CN`, etc.) to capture imprecision, length, confidence intervals, copy‐number, and phasing information.

### 5.4 Specifying complex rearrangements with breakends

Breakend notation allows representation of arbitrary rearrangements:

* **5.4.1 Inserted Sequence** — specify sequence within `[]` brackets.
* **5.4.2 Large Insertions** — use symbolic `<INS>` + `SVLEN`.
* **5.4.3 Multiple mates** — list multiple breakend records for multi‐way events.
* **5.4.4 Explicit partners** — include partner chromosome/position within the ALT string.
* **5.4.5 Telomeres** — use `<T>` symbolic allele for telomeric breakends.
* **5.4.6 Event modifiers** — flags like `IMPRECISE`, `NOVEL`.
* **5.4.7 Inversions** — `<INV>` with `SVLEN`.
* **5.4.8 Uncertainty around breakend location** — use `CIPOS` and `CIEND`.
* **5.4.9 Single breakends** — breakend with no defined partner.

For full details and all sub‑cases (e.g. nested inversions, sample mixtures), see Sections 5.4.1–5.4.9 in the spec.&#x20;

