"""
GFF/GTF Genome Feature Loader for building genomic structures from GFF/GTF files.
"""
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

from .base_parser import BaseParser
from .gtf_parser import GtfParser
from .gff3_parser import Gff3Parser
from .gff2_parser import Gff2Parser
from .format import Format
from .record import GffRecord
from ...genomics.genome import Genome
from ...genomics.chromosome import Chromosome
from ...genomics.gene import Gene
from ...genomics.transcript import Transcript
from ...genomics.exon import Exon
from ...genomics.cds import CDS
from ...genomics.biotype import Biotype
from ...genomics.phase import Phase


class GffGenomeLoader:
    """
    Class for loading genomic features from GFF/GTF files and building a genome structure.
    
    This class parses GFF/GTF files and builds a hierarchical structure of genomic features
    (genes, transcripts, exons, CDS) that can be added to a Genome object.
    """
    
    def __init__(self):
        """Initialize loader with optional progress bars."""
        self.genes_by_id: dict[str, Gene] = {}
        self.transcripts_by_id: dict[str, Transcript] = {}
        self.exons_by_transcript: dict[str, list[Exon]] = defaultdict(list)
        self.cds_by_transcript: dict[str, list[CDS]] = defaultdict(list)
    
    def load(self, gff_file: Path, verbose: bool = None) -> dict[str, Gene]:
        """Load genomic features from a GFF/GTF file."""
        # Use provided verbose parameter if given, otherwise use instance variable
        verbose = self.verbose if verbose is None else verbose
        
        # Reset collections
        self.genes_by_id = {}
        self.transcripts_by_id = {}
        self.exons_by_transcript = defaultdict(list)
        self.cds_by_transcript = defaultdict(list)
        
        # Detect format and get appropriate parser
        parser = self._get_parser_for_file(gff_file)
        
        # First pass: collect all records
        if verbose:
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
        self._build_hierarchy(verbose)
        
        return self.genes_by_id
    
    def _get_parser_for_file(self, gff_file: Path) -> BaseParser:
        """Get the appropriate parser for the file format."""
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
        """Process a single GFF/GTF record."""
        feature_type = record.type
        
        if feature_type == "gene":
            self._process_gene_record(record)
        elif feature_type == "transcript" or feature_type == "mRNA":
            self._process_transcript_record(record)
        elif feature_type == "exon":
            self._process_exon_record(record)
        elif feature_type == "CDS":
            self._process_cds_record(record)
    
    def _process_gene_record(self, record: GffRecord) -> None:
        """Process a gene record."""
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
        """Process a transcript record."""
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
            utrs=[],
            introns=[],
            splice_sites=[]
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
        """Process an exon record."""
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
        """Process a CDS record."""
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
    
    def _build_hierarchy(self, verbose: bool = False) -> None:
        """Build the hierarchical structure of genomic features."""
        # Add exons and CDS to transcripts
        for transcript_id, transcript in self.transcripts_by_id.items():
            transcript.exons = self.exons_by_transcript.get(transcript_id, [])
            transcript.cds_list = self.cds_by_transcript.get(transcript_id, [])
        
        # Add transcripts to genes
        for gene_id, gene in self.genes_by_id.items():
            gene.transcripts = [t for t in self.transcripts_by_id.values() if t.gene_id == gene_id]
    