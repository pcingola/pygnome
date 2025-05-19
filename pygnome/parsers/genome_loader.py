"""
Genome Loader for building complete genomes from annotation and sequence files.
"""
from pathlib import Path
from collections import defaultdict
from enum import Enum
from tqdm import tqdm

from ..genomics.genome import Genome
from ..genomics.chromosome import Chromosome
from ..genomics.gene import Gene
from ..genomics.transcript import Transcript
from ..genomics.exon import Exon
from ..genomics.cds import CDS
from ..genomics.biotype import Biotype
from ..genomics.phase import Phase
from ..genomics.utr import UTR, UTRType
from ..sequences.dna_string import DnaString
from .fasta.fasta_parser import FastaParser
from .gff.base_parser import BaseParser
from .gff.gtf_parser import GtfParser
from .gff.gff3_parser import Gff3Parser
from .gff.gff2_parser import Gff2Parser
from .gff.format import Format
from .gff.record import GffRecord


class ErrorHandling(str, Enum):
    """Enum for error handling options."""
    THROW = "throw"
    WARN = "warn"
    IGNORE = "ignore"


class GenomeLoader:
    """
    Class for loading complete genomes from annotation and sequence files.
    
    This class combines sequence data from FASTA files with annotation data
    from GFF/GTF files to build a complete Genome object with chromosomes,
    genes, transcripts, exons, and other genomic features.
    
    A GenomeLoader object is used to load a single genome, with the GTF file
    and FASTA file provided in the constructor.
    """
    
    def __init__(
        self, 
        annotation_file: Path,
        sequence_file: Path,
        genome_name: str = "genome", 
        species: str | None = None, 
        verbose: bool = False,
        error_handling: ErrorHandling = ErrorHandling.WARN
    ):
        """
        Initialize the genome loader.
        
        Args:
            annotation_file: Path to the GFF/GTF annotation file
            sequence_file: Path to the FASTA sequence file
            genome_name: Name of the genome
            species: Species name
            verbose: Whether to print progress information during loading
            error_handling: How to handle consistency errors (throw, warn, or ignore)
        """
        self.annotation_file = annotation_file
        self.sequence_file = sequence_file
        self.genome_name = genome_name
        self.species = species
        self.verbose = verbose
        self.error_handling = error_handling
        
        # Initialize collections for genomic features
        self.genes_by_id: dict[str, Gene] = {}
        self.transcripts_by_id: dict[str, Transcript] = {}
        self.exons_by_transcript: dict[str, list[Exon]] = defaultdict(list)
        self.cds_by_transcript: dict[str, list[CDS]] = defaultdict(list)
        self.utrs_by_transcript: dict[str, list[UTR]] = defaultdict(list)
        
        # Initialize the genome
        self.genome = Genome(name=genome_name, species=species)
        
    def load(self) -> Genome:
        """
        Load a genome from the annotation and sequence files.
        
        Returns:
            Loaded Genome object
        """
        # Load chromosome sequences
        if self.verbose:
            print(f"Loading chromosome sequences from {self.sequence_file}")
        chromosomes = self.load_sequences(self.sequence_file)
        
        # Load genomic features
        if self.verbose:
            print(f"Loading genomic features from {self.annotation_file}")
        self.load_features(self.annotation_file)
        
        # Add features to chromosomes and add chromosomes to genome
        for chrom in chromosomes.values():
            chrom_genes = [gene for gene in self.genes_by_id.values() if gene.chrom == chrom.name]
            chrom.add_genes(chrom_genes)
            self.genome.add_chromosome(chrom)
            
        # Check for consistency errors if requested
        if self.error_handling != ErrorHandling.IGNORE:
            self.check_consistency()
            
        if self.verbose:
            print(f"Genome loading complete: {len(self.genome.chromosomes)} chromosomes, {len(self.genome.genes)} genes")
        
        return self.genome

    def load_sequences(self, sequence_file: Path) -> dict[str, Chromosome]:
        """
        Load chromosome sequences from a FASTA file.
        
        Args:
            sequence_file: Path to the FASTA file
            
        Returns:
            Dictionary of chromosome names to Chromosome objects
        """
        # Load sequences using the FASTA parser
        parser = FastaParser(sequence_file)
        records = parser.load()
        
        if not records:
            raise ValueError(f"No sequences found in {sequence_file}")
        
        # Create chromosomes from sequences
        chromosomes = {}
        
        if self.verbose and len(records) > 1:
            # Use progress bar for multiple chromosomes
            chrom_iter = tqdm(records, desc="Loading chromosomes")
        else:
            chrom_iter = records
            
        for record in chrom_iter:
            chrom_name = record.identifier
            sequence = record.sequence
            
            # Convert to DnaString if it's a string
            if isinstance(sequence, str):
                sequence = DnaString(sequence)
                
            chromosome = Chromosome(name=chrom_name, sequence=sequence)
            chromosomes[chrom_name] = chromosome
        
        return chromosomes
    
    def load_features(self, gff_file: Path) -> None:
        """
        Load genomic features from a GFF/GTF file.
        
        Args:
            gff_file: Path to the GFF/GTF file
        """
        # Reset collections
        self.genes_by_id = {}
        self.transcripts_by_id = {}
        self.exons_by_transcript = defaultdict(list)
        self.cds_by_transcript = defaultdict(list)
        self.utrs_by_transcript = defaultdict(list)
        
        # Detect format and get appropriate parser
        parser = self._get_parser_for_file(gff_file)
        
        # First pass: collect all records
        if self.verbose:
            # Count lines for progress bar
            line_count = sum(1 for line in open(gff_file, 'r')
                           if line.strip() and not line.startswith('#'))
            
            # Parse with progress bar
            with tqdm(total=line_count, desc="Parsing features") as pbar:
                for record in parser.parse(str(gff_file)):
                    self._process_record(record)
                    pbar.update(1)
        else:
            # Parse without progress bar
            for record in parser.parse(str(gff_file)):
                self._process_record(record)
        
        # Second pass: build the hierarchy
        self._build_hierarchy()
    
    def _get_parser_for_file(self, gff_file: Path) -> BaseParser:
        """
        Get the appropriate parser for the file format.
        
        Args:
            gff_file: Path to the GFF/GTF file
            
        Returns:
            Appropriate parser for the file format
        """
        # Detect format
        format_type = BaseParser.detect_format(str(gff_file))
        
        # Create appropriate parser
        if format_type == Format.GTF:
            return GtfParser()
        elif format_type == Format.GFF3:
            return Gff3Parser()
        else:  # Format.GFF2
            return Gff2Parser()
    
    def _process_record(self, record: GffRecord) -> None:
        """
        Process a single GFF/GTF record.
        
        Args:
            record: GFF/GTF record to process
        """
        feature_type = record.type
        
        if feature_type == "gene":
            self._process_gene_record(record)
        elif feature_type == "transcript" or feature_type == "mRNA":
            self._process_transcript_record(record)
        elif feature_type == "exon":
            self._process_exon_record(record)
        elif feature_type == "CDS":
            self._process_cds_record(record)
        elif feature_type == "five_prime_utr" or feature_type == "5UTR":
            self._process_utr_record(record, UTRType.FIVE_PRIME)
        elif feature_type == "three_prime_utr" or feature_type == "3UTR":
            self._process_utr_record(record, UTRType.THREE_PRIME)
    
    def _process_gene_record(self, record: GffRecord) -> None:
        """
        Process a gene record.
        
        Args:
            record: Gene record to process
        """
        # Handle different attribute names in different formats
        gene_id = (
            record.get_attribute("gene_id") or 
            record.get_attribute("ID") or
            f"gene_{record.chrom}_{record.start}_{record.end}"
        )
        
        gene_name = (
            record.get_attribute("gene_name") or 
            record.get_attribute("Name") or 
            gene_id
        )
        
        gene_biotype = (
            record.get_attribute("gene_type") or 
            record.get_attribute("gene_biotype") or 
            record.get_attribute("biotype") or 
            "unknown"
        )
        
        # Create a Gene object
        gene = Gene(
            id=gene_id,
            chrom=record.chrom,
            start=record.start - 1,  # Convert to 0-based
            end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
            strand=record.strand,
            name=gene_name,
            biotype=Biotype(gene_biotype) if gene_biotype else None,
            transcripts=[]
        )
        self.genes_by_id[gene_id] = gene
    
    def _process_transcript_record(self, record: GffRecord) -> None:
        """
        Process a transcript record.
        
        Args:
            record: Transcript record to process
        """
        # Handle different attribute names in different formats
        gene_id = (
            record.get_attribute("gene_id") or 
            record.get_attribute("Parent") or
            f"gene_{record.chrom}_{record.start}_{record.end}"
        )
        
        transcript_id = (
            record.get_attribute("transcript_id") or 
            record.get_attribute("ID") or
            f"transcript_{record.chrom}_{record.start}_{record.end}"
        )
        
        transcript_biotype = (
            record.get_attribute("transcript_type") or 
            record.get_attribute("transcript_biotype") or 
            record.get_attribute("biotype") or 
            "unknown"
        )
        
        # Create a Transcript object
        transcript = Transcript(
            id=transcript_id,
            chrom=record.chrom,
            start=record.start - 1,  # Convert to 0-based
            end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
            strand=record.strand,
            gene_id=gene_id,
            biotype=Biotype(transcript_biotype) if transcript_biotype else None,
            exons=[],
            cds_list=[],
            utrs=[]
        )
        self.transcripts_by_id[transcript_id] = transcript
        
        # If we don't have a gene record for this transcript, create one
        if gene_id not in self.genes_by_id:
            gene = Gene(
                id=gene_id,
                chrom=record.chrom,
                start=record.start - 1,  # Convert to 0-based
                end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
                strand=record.strand,
                name=gene_id,
                biotype=None,
                transcripts=[]
            )
            self.genes_by_id[gene_id] = gene
    
    def _process_exon_record(self, record: GffRecord) -> None:
        """
        Process an exon record.
        
        Args:
            record: Exon record to process
        """
        # Handle different attribute names in different formats
        transcript_id = (
            record.get_attribute("transcript_id") or 
            record.get_attribute("Parent") or 
            None
        )
        
        if not transcript_id:
            # Skip exons without a parent transcript
            return
        
        exon_id = (
            record.get_attribute("exon_id") or 
            record.get_attribute("ID") or
            f"exon_{record.chrom}_{record.start}_{record.end}"
        )
        
        # Create an Exon object
        exon = Exon(
            id=exon_id,
            chrom=record.chrom,
            start=record.start - 1,  # Convert to 0-based
            end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
            strand=record.strand,
            phase=Phase(record.phase) if record.phase is not None else None
        )
        
        # Handle multiple parent transcripts (GFF3)
        if "," in transcript_id:
            for tid in transcript_id.split(","):
                self.exons_by_transcript[tid].append(exon)
        else:
            self.exons_by_transcript[transcript_id].append(exon)
    
    def _process_cds_record(self, record: GffRecord) -> None:
        """
        Process a CDS record.
        
        Args:
            record: CDS record to process
        """
        # Handle different attribute names in different formats
        transcript_id = (
            record.get_attribute("transcript_id") or 
            record.get_attribute("Parent") or 
            None
        )
        
        if not transcript_id:
            # Skip CDS without a parent transcript
            return
        
        # Create a CDS object
        cds = CDS(
            id=f"CDS_{record.chrom}_{record.start}_{record.end}",
            chrom=record.chrom,
            start=record.start - 1,  # Convert to 0-based
            end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
            strand=record.strand,
            phase=Phase(record.phase) if record.phase is not None else Phase.ZERO
        )
        
        # Handle multiple parent transcripts (GFF3)
        if "," in transcript_id:
            for tid in transcript_id.split(","):
                self.cds_by_transcript[tid].append(cds)
        else:
            self.cds_by_transcript[transcript_id].append(cds)
    
    def _process_utr_record(self, record: GffRecord, utr_type: UTRType) -> None:
        """
        Process a UTR record.
        
        Args:
            record: UTR record to process
            utr_type: Type of UTR (5' or 3')
        """
        # Handle different attribute names in different formats
        transcript_id = (
            record.get_attribute("transcript_id") or 
            record.get_attribute("Parent") or 
            None
        )
        
        if not transcript_id:
            # Skip UTRs without a parent transcript
            return
        
        # Create a UTR object
        utr = UTR(
            id=f"UTR_{utr_type.value}_{record.chrom}_{record.start}_{record.end}",
            chrom=record.chrom,
            start=record.start - 1,  # Convert to 0-based
            end=record.end,  # GFF/GTF is 1-based, inclusive; we want 0-based, exclusive
            strand=record.strand,
            utr_type=utr_type
        )
        
        # Handle multiple parent transcripts (GFF3)
        if "," in transcript_id:
            for tid in transcript_id.split(","):
                self.utrs_by_transcript[tid].append(utr)
        else:
            self.utrs_by_transcript[transcript_id].append(utr)
    
    def _build_hierarchy(self) -> None:
        """Build the hierarchical structure of genomic features."""
        # Add exons, CDS, and UTRs to transcripts
        for transcript_id, transcript in self.transcripts_by_id.items():
            transcript.exons = self.exons_by_transcript.get(transcript_id, [])
            transcript.cds_list = self.cds_by_transcript.get(transcript_id, [])
            transcript.utrs = self.utrs_by_transcript.get(transcript_id, [])
        
        # Add transcripts to genes
        for gene_id, gene in self.genes_by_id.items():
            gene.transcripts = [t for t in self.transcripts_by_id.values() if t.gene_id == gene_id]
    
    def check_consistency(self) -> None:
        """
        Check for consistency errors in the loaded genome.
        
        This method checks for various consistency issues, such as:
        - Transcripts without exons
        - CDS regions outside of exons
        - UTRs outside of exons
        - Overlapping exons
        - Overlapping CDS regions
        
        Depending on the error_handling setting, issues will either:
        - Raise an exception (THROW)
        - Print a warning message (WARN)
        - Be silently ignored (IGNORE)
        """
        errors = []
        
        # Check each transcript
        for transcript_id, transcript in self.transcripts_by_id.items():
            # Check if transcript has exons
            if not transcript.exons:
                errors.append(f"Transcript {transcript_id} has no exons")
                continue
            
            # Check if CDS regions are within exons
            for cds in transcript.cds_list:
                if not any(exon.start <= cds.start and exon.end >= cds.end for exon in transcript.exons):
                    errors.append(f"CDS {cds.id} in transcript {transcript_id} is not within any exon")
            
            # Check if UTRs are within exons
            for utr in transcript.utrs:
                if not any(exon.start <= utr.start and exon.end >= utr.end for exon in transcript.exons):
                    errors.append(f"UTR {utr.id} in transcript {transcript_id} is not within any exon")
            
            # Check for overlapping exons
            sorted_exons = sorted(transcript.exons, key=lambda x: x.start)
            for i in range(len(sorted_exons) - 1):
                if sorted_exons[i].end > sorted_exons[i + 1].start:
                    errors.append(f"Overlapping exons in transcript {transcript_id}: {sorted_exons[i].id} and {sorted_exons[i + 1].id}")
            
            # Check for overlapping CDS regions
            sorted_cds = sorted(transcript.cds_list, key=lambda x: x.start)
            for i in range(len(sorted_cds) - 1):
                if sorted_cds[i].end > sorted_cds[i + 1].start:
                    errors.append(f"Overlapping CDS regions in transcript {transcript_id}: {sorted_cds[i].id} and {sorted_cds[i + 1].id}")
        
        # Handle errors based on error_handling setting
        if errors and self.error_handling == ErrorHandling.THROW:
            raise ValueError("\n".join(errors))
        elif errors and self.error_handling == ErrorHandling.WARN:
            for error in errors:
                print(f"WARNING: {error}")
    
    def __str__(self) -> str:
        """
        Return a string representation of the genome loader.
        
        This method provides a summary of the loaded genome without
        attempting to print the entire genome data, which could be
        several gigabytes in size.
        """
        genome_str = f"GenomeLoader({self.genome_name}"
        if self.species:
            genome_str += f", {self.species}"
        genome_str += ")"
        
        if not self.genome.chromosomes:
            return f"{genome_str} - No chromosomes loaded"
        
        chrom_count = len(self.genome.chromosomes)
        gene_count = len(self.genome.genes)
        
        # Calculate total sequence length
        total_length = 0
        for chrom in self.genome.chromosomes.values():
            if chrom.length is not None:
                total_length += chrom.length
        
        # Format total length in a human-readable way
        if total_length > 1_000_000_000:
            length_str = f"{total_length / 1_000_000_000:.2f} Gb"
        elif total_length > 1_000_000:
            length_str = f"{total_length / 1_000_000:.2f} Mb"
        elif total_length > 1_000:
            length_str = f"{total_length / 1_000:.2f} kb"
        else:
            length_str = f"{total_length} bp"
        
        return f"{genome_str} - {chrom_count} chromosomes, {gene_count} genes, {length_str}"