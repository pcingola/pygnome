# GFF2

Reference: https://gmod.org/wiki/GFF2

[GFF2](http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml) is a supported format in
GMOD, **but it is now deprecated and if you have a choice you should use**
**[GFF3](https://gmod.org/wiki/GFF3 "GFF3")**. Unfortunately, data is sometimes only available
in GFF2 format. GFF2 has a number of shortcomings compared to GFF3. GFF2
can only represent 2 level feature hierarchies, while GFF3 can support
arbitrary levels. GFF2 also does not require that column 3, the feature
type, be part of the sequence ontology. It can be any string. This often
led to quality control and data exchange problems.

## Contents [Anchor](https://gmod.org/wiki/GFF2\#contents)

- [1GFF2 is\\
Deprecated!](https://gmod.org/wiki/GFF2#GFF2_is_Deprecated.21)
  - [1.1Why GFF2\\
    is harmful to your\\
    health](https://gmod.org/wiki/GFF2#Why_GFF2_is_harmful_to_your_health)
- [2The GFF2 File\\
Format](https://gmod.org/wiki/GFF2#The_GFF2_File_Format)
  - [2.1Creating a\\
    GFF2 table](https://gmod.org/wiki/GFF2#Creating_a_GFF2_table)
    - [2.1.1Using\\
      the Group field for simple\\
      features](https://gmod.org/wiki/GFF2#Using_the_Group_field_for_simple_features)
    - [2.1.2Using\\
      the Group field to group features that belong\\
      together](https://gmod.org/wiki/GFF2#Using_the_Group_field_to_group_features_that_belong_together)
    - [2.1.3Using\\
      the Group field to add a\\
      note](https://gmod.org/wiki/GFF2#Using_the_Group_field_to_add_a_note)
    - [2.1.4Using\\
      the Group field to add an alternative\\
      name](https://gmod.org/wiki/GFF2#Using_the_Group_field_to_add_an_alternative_name)
  - [2.2Identifying the reference\\
    sequence](https://gmod.org/wiki/GFF2#Identifying_the_reference_sequence)
  - [2.3Sequence\\
    alignments](https://gmod.org/wiki/GFF2#Sequence_alignments)
  - [2.4Dense\\
    quantitative data](https://gmod.org/wiki/GFF2#Dense_quantitative_data)
  - [2.5Loading\\
    the GFF file into the\\
    database](https://gmod.org/wiki/GFF2#Loading_the_GFF_file_into_the_database)
  - [2.6Aggregators](https://gmod.org/wiki/GFF2#Aggregators)
- [3Converting\\
GFF2 to GFF3](https://gmod.org/wiki/GFF2#Converting_GFF2_to_GFF3)
  - [3.1Column 3:\\
    Feature Type](https://gmod.org/wiki/GFF2#Column_3:_Feature_Type)
  - [3.2Column 9:\\
    Group / Attributes](https://gmod.org/wiki/GFF2#Column_9:_Group_.2F_Attributes)

## GFF2 is Deprecated! [Anchor](https://gmod.org/wiki/GFF2\#gff2-is-deprecated)

The GFF file format stands for “Gene Finding Format” or or “General
Feature Format” and was invented at the Sanger Centre. It is easy to
use, but it suffers from two main limitations (see the box).

### Why GFF2 is harmful to your health [Anchor](https://gmod.org/wiki/GFF2\#why-gff2-is-harmful-to-your-health)

One of GFF2’s problems is that it is only able to represent one level of
nesting of features. This is mainly a problem when dealing with genes
that have multiple alternatively-spliced transcripts. GFF2 is unable to
deal with the three-level hierarchy of _gene → transcript → exon_. Most
people get around this by declaring a series of transcripts and giving
them similar names to indicate that they come from the same gene. The
second limitation is that while GFF2 allows you to create two-level
hierarchies, such as _transcript → exon_, it doesn’t have any concept of
the direction of the hierarchy. So it doesn’t know whether the exon is a
subfeature of the transcript, or vice-versa. This means you have to use
“aggregators” to sort out the relationships. This is a major pain in the
neck. For this reason, GFF2 format has been deprecated in favor of
[GFF3](https://gmod.org/wiki/GFF3 "GFF3") format databases.

See [GFF3](https://gmod.org/wiki/GFF3 "GFF3") for more on the current version of GFF.

## The GFF2 File Format [Anchor](https://gmod.org/wiki/GFF2\#the-gff2-file-format)

The GFF format is a flat tab-delimited file, each line of which
corresponds to an annotation, or feature. Each line has nine columns and
looks like this:

```
Chr1  curated  CDS 365647  365963  .  +  1  Transcript "R119.7"

```

The 9 columns are as follows:

reference sequence

This is the ID of the sequence that is used to establish the coordinate
system of the annotation. In the example above, the reference sequence
is “Chr1”.

source

The source of the annotation. This field describes how the annotation
was derived. In the example above, the source is “curated” to indicate
that the feature is the result of human curation. The names and versions
of software programs are often used for the source field, as in
“tRNAScan-SE/1.2”.

method

The annotation method, also known as type. This field describes the type
of the annotation, such as “CDS”. Together the method and source
describe the annotation type.

start position

The start of the annotation relative to the reference sequence.

stop position

The stop of the annotation relative to the reference sequence. Start is
always less than or equal to stop.

score

For annotations that are associated with a numeric score (for example, a
sequence similarity), this field describes the score. The score units
are completely unspecified, but for sequence similarities, it is
typically percent identity. Annotations that do not have a score can use
“.”

strand

For those annotations which are strand-specific, this field is the
strand on which the annotation resides. It is “+” for the forward
strand, “-“ for the reverse strand, or “.” for annotations that are not
stranded.

phase

For annotations that are linked to proteins, this field describes the
phase of the annotation on the codons. It is a number from 0 to 2, or
“.” for features that have no phase.

group

GFF provides a simple way of generating annotation hierarchies (“is
composed of” relationships) by providing a group field. The group field
contains the class and ID of an annotation which is the logical parent
of the current one. In the example given above, the group is the
Transcript named “R119.7”.

The group field is also used to store information about the target of
sequence similarity hits, and miscellaneous notes. See the next section
for a description of how to describe similarity targets.

The sequences used to establish the coordinate system for annotations
can correspond to sequenced clones, clone fragments, contigs or
super-contigs.

In addition to a group ID, the GFF format allows annotations to have a
group class. This makes sure that all groups are unique even if they
happen to share the same name. For example, you can have a GenBank
accession named AP001234 and a clone named AP001234 and distinguish
between them by giving the first one a class of Accession and the second
a class of Clone.

You should use double-quotes around the group name or class if it
contains white space.

### Creating a GFF2 table [Anchor](https://gmod.org/wiki/GFF2\#creating-a-gff2-table)

The first 8 fields of the GFF2 format are easy to understand. The group
field is a challenge. It is used in several distinct ways:

- to group together a single sequence feature that spans a discontinuous
range, such as a gapped alignment.
- to name a feature, allowing it to be retrieved by name.
- to add one or more notes to the annotation.
- to add an alternative name

#### Using the Group field for simple features [Anchor](https://gmod.org/wiki/GFF2\#using-the-group-field-for-simple-features)

For a simple feature that spans a single continuous range, choose a name
and class for the object and give it a line in the GFF2 file that refers
to its start and stop positions.

```
Chr3   giemsa heterochromatin  4500000 6000000 . . .   Band 3q12.1

```

#### Using the Group field to group features that belong together [Anchor](https://gmod.org/wiki/GFF2\#using-the-group-field-to-group-features-that-belong-together)

For a group of features that belong together, such as the exons in a
transcript, choose a name and class for the object. Give each segment a
separate line in the GFF2 file but use the same name for each line. For
example:

```
IV     curated exon    5506900 5506996 . + .   Transcript B0273.1
IV     curated exon    5506026 5506382 . + .   Transcript B0273.1
IV     curated exon    5506558 5506660 . + .   Transcript B0273.1
IV     curated exon    5506738 5506852 . + .   Transcript B0273.1

```

These four lines refer to a biological object of class “Transcript” and
name B0273.1. Each of its parts uses the method “exon”, source
“curated”. Once loaded, the user will be able to search the genome for
this object by asking the browser to retrieve “Transcript:B0273.1”. The
browser can also be configured to allow the Transcript: prefix to be
omitted.

You can extend the idiom for objects that have heterogeneous parts, such
as a transcript that has 5’ and 3’ UTRs

```
IV     curated  mRNA   5506800 5508917 . + .   Transcript B0273.1; Note "Zn-Finger"
IV     curated  5'UTR  5506800 5508999 . + .   Transcript B0273.1
IV     curated  exon   5506900 5506996 . + .   Transcript B0273.1
IV     curated  exon   5506026 5506382 . + .   Transcript B0273.1
IV     curated  exon   5506558 5506660 . + .   Transcript B0273.1
IV     curated  exon   5506738 5506852 . + .   Transcript B0273.1
IV     curated  3'UTR  5506852 5508917 . + .   Transcript B0273.1

```

In this example, there is a single feature with method “mRNA” that spans
the entire range. It is grouped with subparts of type 5’UTR, 3’UTR and
exon. They are all grouped together into a Transcript named B0273.1.
Furthermore the mRNA feature has a note attached to it.

**NOTE:** The subparts of a feature are in absolute (chromosomal or
contig) coordinates. It is not currently possible to define a feature in
absolute coordinates and then to load its subparts using coordinates
that are relative to the start of the feature.

Some annotations do not need to be individually named. For example, it
is probably not useful to assign a unique name to each ALU repeat in a
vertebrate genome. For these, just leave the Group field empty.

#### Using the Group field to add a note [Anchor](https://gmod.org/wiki/GFF2\#using-the-group-field-to-add-a-note)

The group field can be used to add one or more notes to an annotation.
To do this, place a semicolon after the group name and add a Note field:

```
Chr3 giemsa heterochromatin 4500000 6000000 . . . Band 3q12.1 ; Note "Marfan's syndrome"

```

You can add multiple Notes. Just separate them by semicolons:

```
 Band 3q12.1 ; Note "Marfan's syndrome" ; Note "dystrophic dysplasia"

```

The Note should come AFTER the group type and name.

#### Using the Group field to add an alternative name [Anchor](https://gmod.org/wiki/GFF2\#using-the-group-field-to-add-an-alternative-name)

If you want the feature to be quickly searchable by an alternative name,
you can add one or more Alias tags. A feature can have multiple aliases,
and multiple features can share the same alias:

```
Chr3 giemsa heterochromatin 4500000 6000000 . . . Band 3q12.1 ; Alias MFX

```

Searches for aliases will be both faster and more reliable than searches
for keywords in notes, since the latter relies on whole-text search
methods that vary somewhat from DBMS to DBMS.

### Identifying the reference sequence [Anchor](https://gmod.org/wiki/GFF2\#identifying-the-reference-sequence)

Each reference sequence in the GFF table must itself have an entry. This
is necessary so that the length of the reference sequence is known.

For example, if “Chr1” is used as a reference sequence, then the GFF
file should have an entry for it similar to this one:

```
Chr1 assembly chromosome 1 14972282 . + . Sequence Chr1

```

This indicates that the reference sequence named “Chr1” has length
14972282 bp, method “chromosome” and source “assembly”. In addition, as
indicated by the group field, Chr1 has class “Sequence” and name “Chr1”.

It is suggested that you use “Sequence” as the class name for all
reference sequences, since this is the default class used by the
Bio::DB::GFF module when no more specific class is requested. If you use
a different class name, then be sure to indicate that fact with the
“reference class” option (see below).

### Sequence alignments [Anchor](https://gmod.org/wiki/GFF2\#sequence-alignments)

There are several cases in which an annotation indicates the
relationship between two sequences. One common one is a similarity hit,
where the annotation indicates an alignment. A second common case is a
map assembly, in which the annotation indicates that a portion of a
larger sequence is built up from one or more smaller ones.

Both cases are indicated by using the Target tag in the group field. For
example, a typical similarity hit will look like this:

```
Chr1 BLASTX similarity 76953 77108 132 + 0 Target Protein:SW:ABL_DROME 493 544

```

Here, the group field contains the Target tag, followed by an identifier
for the biological object. The GFF format uses the notation Class:Name
for the biological object, and even though this is stylistically
inconsistent, that’s the way it’s done. The object identifier is
followed by two integers indicating the start and stop of the alignment
on the target sequence.

Unlike the main start and stop columns, it is possible for the target
start to be greater than the target end. The previous example indicates
that the the section of Chr1 from 76,953 to 77,108 aligns to the protein
SW:ABL\_DROME starting at position 493 and extending to position 544.

A similar notation is used for sequence assembly information as shown in
this example:

```
Chr1        assembly Link   10922906 11177731 . . . Target Sequence:LINK_H06O01 1 254826
LINK_H06O01 assembly Cosmid 32386    64122    . . . Target Sequence:F49B2       6 31742

```

This indicates that the region between bases 10922906 and 11177731 of
Chr1 are composed of LINK\_H06O01 from bp 1 to bp 254826. The region of
LINK\_H0601 between 32386 and 64122 is, in turn, composed of the bases 5
to 31742 of cosmid F49B2.

### Dense quantitative data [Anchor](https://gmod.org/wiki/GFF2\#dense-quantitative-data)

If you have dense quantitative data, such as tiling array data,
microarray expression data, ChIP-chip or ChIP-seq chromatin
immunoprecipitation data, then you will probably want to create “Wiggle”
format binary files, which represent the quantitative data in a compact
format in external files. Use the `wiggle2gff3.pl` script, included in
this distribution, to format and load this data. Run `wiggle2gff3.pl -h`
for instructions.

### Loading the GFF file into the database [Anchor](https://gmod.org/wiki/GFF2\#loading-the-gff-file-into-the-database)

Use the [BioPerl](https://gmod.org/wiki/BioPerl "BioPerl") script utilities
`bp_bulk_load_gff.pl`, `bp_load_gff.pl` or (if you are brave)
`bp_fast_load_gff.pl` to load the GFF file into the database. For
example, if your database is a MySQL database on the local host named
“dicty”, you can load it into an empty database using
`bp_bulk_load_gff.pl` like this:

```
 bp_bulk_load_gff.pl -c -d dicty my_data.gff

```

To update existing databases, use either `bp_load_gff.pl` or
`bp_fast_load_gff.pl`. The latter is somewhat experimental, so use with
care.

### Aggregators [Anchor](https://gmod.org/wiki/GFF2\#aggregators)

It is not necessary to use aggregators with the
[Chado](https://gmod.org/wiki/Chado "Chado"),
[BioSQL](https://gmod.org/wiki/BioSQL "BioSQL"), or Bio::DB::SeqFeature::Store [GBrowse\\
Adaptors](https://gmod.org/wiki/GBrowse_Adaptors "GBrowse Adaptors"), or any other adaptor
that is based on [GFF3](https://gmod.org/wiki/GFF3 "GFF3").

The Bio::DB::GFF adaptor (and only Bio::DB::GFF!) has a feature known as
“aggregators”. These are small software packages that recognize certain
common feature types and convert them into complex biological objects.
These aggregators make it possible to develop intelligent graphical
representations of annotations, such as a gene that draws confirmed
exons differently from predicted ones.

An aggregator typically creates a new composite feature with a different
method than any of its components. For example, the standard “alignment”
aggregator takes multiple alignments of method “similarity”, groups them
by their name, and returns a single feature of method “alignment”.

The various aggregators are described in detail in the Bio::DB::GFF
perldoc page. It is easy to write new aggregators, and also possible to
define aggregators on the fly in the GBrowse configuration file. It is
suggested that you use the sample GFF2 files from the yeast,
_drosophila_ and _C. elegans_ projects to see what methods to use to
achieve the desired results.

In addition to the standard aggregators that are distributed with
[BioPerl](https://gmod.org/wiki/BioPerl "BioPerl"), [GBrowse](https://gmod.org/wiki/GBrowse.1 "GBrowse") distributes
several experimental and/or special-purpose aggregators:

match\_gap

This aggregator is used for GFF3 style gapped alignments, in which there
is a single feature of method ‘match’ with a ‘Gap’ attribute. This
aggregator was contributed by Dmitri Bichko.

orf

This aggregator aggregates raw “ORF” features into “coding” features. It
is basically identical to the “coding” aggregator, except that it looks
for features of type “ORF” rather than “cds”.

reftranscript

This aggregator was written to make the compound feature,
“reftranscript” for use with GBrowse editing software developed outside
of the GMOD development group. It can be used to aggregate
“reftranscripts” from “refexons”, loaded as second copy features. These
features, in contrast to “transcripts”, are usually implemented as
features which cannot be edited and serve as starting point references
for annotations added using [GBrowse](https://gmod.org/wiki/GBrowse.1 "GBrowse") for feature
[visualization](https://gmod.org/wiki/Visualization "Visualization"). Adding features to the
compound feature, “reftranscript”, can be done by adding to the
“part\_names” call (i.e. “refCDS”).

waba\_alignment

This aggregator handles the type of alignments produced by Jim Kent’s
WABA program, and was written to be compatible with the _C. elegans_
GFF2 files. It aggregates the following feature types into an aggregate
type of “waba\_alignment”:

- nucleotide\_match:waba\_weak
- nucleotide\_match:waba\_strong
- nucleotide\_match:waba\_coding

wormbase\_gene

This aggregator was written to be compatible with the _C. elegans_ GFF2
files distributed by the Sanger Institute. It aggregates raw “CDS”,
“5’UTR”, “3’UTR”, “polyA” and “TSS” features into “transcript” features.
For compatibility with the idiosyncrasies of the Sanger GFF2 format, it
expects that the full range of the transcript is contained in a main
feature of type “Sequence”.

It is strongly recommended that for mirroring _C. elegans_ annotations,
you use the “processed\_transcript” aggregator in conjunction with the
[GFF3](https://gmod.org/wiki/GFF3 "GFF3") files found at:

[ftp://ftp.wormbase.org/pub/wormbase/genomes/elegans/genome\_feature\_tables/GFF3](ftp://ftp.wormbase.org/pub/wormbase/genomes/elegans/genome_feature_tables/GFF3)

## Converting GFF2 to GFF3 [Anchor](https://gmod.org/wiki/GFF2\#converting-gff2-to-gff3)

Converting a file from GFF2 to [GFF3](https://gmod.org/wiki/GFF3 "GFF3") format is problematic
for several reasons. However, there are several GFF2 to GFF3 converters
available on the web, but each makes specific assumptions about the GFF2
data that limit its applicability. GMOD does not endorse (or disparage)
any particular converter. If you have GFF2 data from an external source,
and they don’t also provide it in GFF3 format, then you may be stuck
with GFF2.

Some areas that need to be addressed by any GFF2 to GFF3 converter:

### Column 3: Feature Type [Anchor](https://gmod.org/wiki/GFF2\#column-3-feature-type)

If the GFF2 file does not use Sequence Ontology terms in column 3 then
some sort of translation will need to be done on the types in the GFF2
to convert them to be SO terms.

### Column 9: Group / Attributes [Anchor](https://gmod.org/wiki/GFF2\#column-9-group--attributes)

Column 9 has a slightly different format and is much more tightly
defined in GFF3 than GFF2. Both require attention. GFF2 does not have
any reserved attribute names, uses C style encoding/escaping of special
characters, and has many other small differences.

Another big problem is that GFF2 supports only one level of feature
nesting. While you can certainly reproduce this minimal nesting in GFF3,
it would be better to also convert your feature representations to be
multi-level at the time you migrate the data to GFF3. This is
non-trivial.

[Categories](https://gmod.org/wiki/Special%253ACategories "Special%253ACategories"):

- [Annotation](https://gmod.org/wiki/Category%253AAnnotation "Category%253AAnnotation")
- [Computing](https://gmod.org/wiki/Category%253AComputing "Category%253AComputing")

## Navigation menu [Anchor](https://gmod.org/wiki/GFF2\#navigation-menu)

### Navigation [Anchor](https://gmod.org/wiki/GFF2\#navigation)

- [GMOD Home](https://gmod.org/wiki/Main_Page)
- [Software](https://gmod.org/wiki/GMOD_Components)
- [Categories /\\
Tags](https://gmod.org/wiki/Categories)

### Documentation [Anchor](https://gmod.org/wiki/GFF2\#documentation)

- [Overview](https://gmod.org/wiki/Overview)
- [FAQs](https://gmod.org/wiki/Category%253AFAQ)
- [HOWTOs](https://gmod.org/wiki/Category%253AHOWTO)
- [Glossary](https://gmod.org/wiki/Glossary)

### Community [Anchor](https://gmod.org/wiki/GFF2\#community)

- [GMOD News](https://gmod.org/wiki/GMOD_News)
- [Training /\\
Outreach](https://gmod.org/wiki/Training_and_Outreach)
- [Support](https://gmod.org/wiki/Support)
- [GMOD Promotion](https://gmod.org/wiki/GMOD_Promotion)
- [Meetings](https://gmod.org/wiki/Meetings)
- [Calendar](https://gmod.org/wiki/Calendar)

### Tools [Anchor](https://gmod.org/wiki/GFF2\#tools)

- [Browse properties](https://gmod.org/wiki/Special%253ABrowse/GFF2)

- Last updated at 13:27 on 21 April
2017.
- Content is available under
[a GNU Free Documentation License](http://www.gnu.org/licenses/fdl-1.3.html) unless otherwise
noted.