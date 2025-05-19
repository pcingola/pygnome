"""Transcript class for genomic annotations."""

from typing import Any
from dataclasses import dataclass, field

from .biotype import Biotype
from .cds import CDS
from .codon_table import CodonTable
from .exon import Exon
from .genomic_feature import GenomicFeature
from .intron import Intron
from .phase import Phase
from .splice import SpliceSite, SpliceSiteType, SpliceRegion
from .strand import Strand
from .utr import UTR, UTRType
from ..sequences.utils import reverse_complement


@dataclass
class Transcript(GenomicFeature):
    """A transcript of a gene, composed of exons, introns, UTRs, and CDS segments."""
    gene_id: str
    biotype: Biotype | None = None
    exons: list[Exon] = field(default_factory=list)
    cds_list: list[CDS] = field(default_factory=list)
    utrs: list[UTR] = field(default_factory=list)
    _introns: list[Intron] | None = None
    splice_sites: list[SpliceSite] = field(default_factory=list)
    splice_regions: list[SpliceRegion] = field(default_factory=list)
    _cds: str | None = None
    _mrna: str | None = None
    _protein: str | None = None
    _aa_to_cds_pos: list[int] | None = None
    _cds_to_aa_pos: list[int] | None = None
    _cds_to_genomic_pos: list[int] | None = None
    _genomic_to_cds_pos: dict[int, int] | None = None
    gene: Any = None  # Reference to Gene
    
    @property
    def coding_length(self) -> int:
        """Return the total length of coding sequence."""
        return sum(cds.length for cds in self.cds_list)
    
    @property
    def exonic_length(self) -> int:
        """Return the total length of exons."""
        return sum(exon.length for exon in self.exons)
    
    @property
    def has_cds(self) -> bool:
        """Return True if the transcript has CDS segments."""
        return len(self.cds_list) > 0
    
    @property
    def is_coding(self) -> bool:
        """
        Return True if the transcript is protein-coding.
        
        This is determined by:
        1. Checking if the biotype is PROTEIN_CODING
        2. If biotype is not set, checking if there are CDS segments
        3. If no CDS segments, checking if there are both 5' and 3' UTRs
        """
        # If biotype is explicitly set to PROTEIN_CODING, it's coding
        if self.biotype == Biotype.PROTEIN_CODING:
            return True
        
        # If it has CDS segments, it's coding
        if self.has_cds:
            return True
        
        # If it has both 5' and 3' UTRs, it's likely coding
        if self.five_prime_utrs and self.three_prime_utrs:
            return True
        
        # Otherwise, assume it's not coding
        return False
    
    @property
    def five_prime_utrs(self) -> list[UTR]:
        """Return all 5' UTRs."""
        return [utr for utr in self.utrs if utr.utr_type == UTRType.FIVE_PRIME]
    
    @property
    def three_prime_utrs(self) -> list[UTR]:
        """Return all 3' UTRs."""
        return [utr for utr in self.utrs if utr.utr_type == UTRType.THREE_PRIME]
    
    def __iter__(self):
        """Iterate over exons sorted by start position."""
        return iter(sorted(self.exons, key=lambda x: x.start))
    
    def cds(self, genome_sequence: dict[str, str]) -> str:
        """Get the coding DNA sequence for this transcript."""
        # Use cached value if available
        if self._cds is not None:
            return self._cds
        
        if not self.is_coding:
            self._cds = ""
            return ""
        
        # Get the chromosome sequence
        if self.chrom not in genome_sequence:
            raise ValueError(f"Chromosome {self.chrom} not found in genome sequence")
        
        chrom_seq = genome_sequence[self.chrom]
        
        # Sort CDS segments by position based on strand
        if self.strand == Strand.POSITIVE:
            # For positive strand, sort by increasing position
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start)
        else:
            # For negative strand, sort by decreasing position
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start, reverse=True)
        
        # Extract the coding sequence for each CDS segment
        coding_seq = ""
        for cds in sorted_cds:
            # Extract the sequence for this CDS segment
            segment_seq = chrom_seq[cds.start:cds.end]
            
            # Reverse complement if on negative strand
            if self.strand == Strand.NEGATIVE:
                segment_seq = reverse_complement(segment_seq)
                
            coding_seq += segment_seq
        
        # Cache the result
        self._cds = coding_seq
        
        return coding_seq
    
    def protein(self, genome_sequence: dict[str, str]) -> str:
        """Get the protein sequence for this transcript."""
        # Get the coding sequence
        coding_seq = self.cds(genome_sequence)
        
        if not coding_seq:
            return ""
        
        # Use the genome's codon table
        table_type = self.gene.chromosome.genome.codon_table_type
        
        # Create a codon table
        codon_table = CodonTable(table_type)
        
        # Translate the coding sequence
        protein_seq = codon_table.translate_sequence(coding_seq)
        
        return protein_seq
    
    def cds_to_aa_pos(self, cds_pos: int) -> int:
        """
        Convert a CDS position to an amino acid position.
        
        Args:
            cds_pos: Position in the CDS (0-based)
            
        Returns:
            Amino acid position (0-based)
        """
        if cds_pos < 0:
            raise ValueError(f"CDS position must be non-negative: {cds_pos}")
        
        # Integer division to get the amino acid position
        return cds_pos // 3
    
    def aa_to_cds_pos(self, aa_pos: int) -> list[int]:
        """
        Convert an amino acid position to CDS positions.
        
        Args:
            aa_pos: Amino acid position (0-based)
            
        Returns:
            List of three CDS positions corresponding to the codon
        """
        if aa_pos < 0:
            raise ValueError(f"Amino acid position must be non-negative: {aa_pos}")
        
        # Calculate the three positions in the codon
        start_pos = aa_pos * 3
        return [start_pos, start_pos + 1, start_pos + 2]
    
    def _build_cds_to_genomic_map(self, genome_sequence: dict[str, str]) -> list[int]:
        """
        Build a mapping from CDS positions to genomic positions.
        
        Args:
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            List where index is CDS position and value is genomic position
        """
        if not self.has_cds:
            return []
        
        # Get the coding sequence to determine length
        coding_seq = self.cds(genome_sequence)
        cds_length = len(coding_seq)
        
        # Initialize the mapping array
        cds_to_genomic = [-1] * cds_length
        
        # Sort CDS segments by position based on strand
        if self.strand == Strand.POSITIVE:
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start)
        else:
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start, reverse=True)
        
        # Fill in the mapping
        cds_pos = 0
        for cds in sorted_cds:
            cds_segment_length = cds.length
            
            if self.strand == Strand.POSITIVE:
                for i in range(cds_segment_length):
                    if cds_pos < cds_length:
                        cds_to_genomic[cds_pos] = cds.start + i
                        cds_pos += 1
            else:
                for i in range(cds_segment_length):
                    if cds_pos < cds_length:
                        cds_to_genomic[cds_pos] = cds.end - 1 - i
                        cds_pos += 1
        
        self._cds_to_genomic_pos = cds_to_genomic
        return cds_to_genomic
    
    def _build_genomic_to_cds_map(self, genome_sequence: dict[str, str]) -> dict[int, int]:
        """
        Build a mapping from genomic positions to CDS positions.
        
        Args:
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            Dictionary mapping genomic positions to CDS positions
        """
        if not self.has_cds:
            return {}
        
        # Get the CDS to genomic mapping
        if self._cds_to_genomic_pos is None:
            self._build_cds_to_genomic_map(genome_sequence)
        
        # Create the reverse mapping
        genomic_to_cds = {}
        
        # Ensure we have a valid mapping
        if self._cds_to_genomic_pos is not None:
            for cds_pos, genomic_pos in enumerate(self._cds_to_genomic_pos):
                if genomic_pos >= 0:
                    genomic_to_cds[genomic_pos] = cds_pos
        
        self._genomic_to_cds_pos = genomic_to_cds
        return genomic_to_cds
    
    def genomic_pos_to_cds_pos(self, genomic_pos: int, genome_sequence: dict[str, str]) -> int:
        """
        Convert a genomic position to a CDS position.
        
        Args:
            genomic_pos: Genomic position
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            CDS position or -1 if the position is not in a CDS
        """
        # Build the mapping if needed
        if self._genomic_to_cds_pos is None:
            self._build_genomic_to_cds_map(genome_sequence)
        
        # Look up the position
        if self._genomic_to_cds_pos is not None:
            return self._genomic_to_cds_pos.get(genomic_pos, -1)
        return -1
    
    def cds_pos_to_genomic_pos(self, cds_pos: int, genome_sequence: dict[str, str]) -> int:
        """
        Convert a CDS position to a genomic position.
        
        Args:
            cds_pos: CDS position
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            Genomic position or -1 if the CDS position is invalid
        """
        # Build the mapping if needed
        if self._cds_to_genomic_pos is None:
            self._build_cds_to_genomic_map(genome_sequence)
        
        # Check if the position is valid
        if self._cds_to_genomic_pos is None:
            return -1
            
        if cds_pos < 0 or cds_pos >= len(self._cds_to_genomic_pos):
            return -1
        
        return self._cds_to_genomic_pos[cds_pos]
    
    def mrna(self, genome_sequence: dict[str, str]) -> str:
        """
        Get the full mRNA sequence (5' UTR + CDS + 3' UTR) for this transcript.
        
        Args:
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            The mRNA sequence as a string
        """
        # Use cached value if available
        if self._mrna is not None:
            return self._mrna
        
        # Get the chromosome sequence
        if self.chrom not in genome_sequence:
            raise ValueError(f"Chromosome {self.chrom} not found in genome sequence")
        
        chrom_seq = genome_sequence[self.chrom]
        
        # Sort exons by position based on strand
        if self.strand == Strand.POSITIVE:
            sorted_exons = sorted(self.exons, key=lambda x: x.start)
        else:
            sorted_exons = sorted(self.exons, key=lambda x: x.start, reverse=True)
        
        # Extract the sequence for each exon
        mrna_seq = ""
        for exon in sorted_exons:
            # Extract the sequence for this exon
            segment_seq = chrom_seq[exon.start:exon.end]
            
            # Reverse complement if on negative strand
            if self.strand == Strand.NEGATIVE:
                segment_seq = reverse_complement(segment_seq)
                
            mrna_seq += segment_seq
        
        # Cache the result
        self._mrna = mrna_seq
        
        return mrna_seq
    
    def infer_phase(self, genome_sequence: dict[str, str]) -> None:
        """
        Infer the phase for CDS segments when missing.
        
        This method calculates the correct phase for each CDS segment based on
        its position in the transcript. It updates the phase attribute of each
        CDS segment in place.
        
        Args:
            genome_sequence: Dictionary mapping chromosome names to sequences
        """
        if not self.has_cds:
            return
        
        # Sort CDS segments by position based on strand
        if self.strand == Strand.POSITIVE:
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start)
        else:
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start, reverse=True)
        
        # The first CDS should have phase 0 (no bases should be skipped)
        current_phase = Phase.ZERO
        
        # Update the phase for each CDS segment
        for i, cds in enumerate(sorted_cds):
            # Set the phase for this CDS segment
            cds.phase = current_phase
            
            # Calculate the phase for the next CDS segment
            # The phase is the number of bases that should be skipped to reach
            # the first base of the next codon
            cds_length = cds.length
            
            # Calculate how many complete codons are in this CDS
            complete_codons = (cds_length - current_phase) // 3
            
            # Calculate how many bases are left over
            bases_left = (cds_length - current_phase) % 3
            
            # The phase for the next CDS is (3 - bases_left) % 3
            # This gives us 0 if bases_left is 0, 1 if bases_left is 2, and 2 if bases_left is 1
            if i < len(sorted_cds) - 1:  # Skip for the last CDS
                next_phase = (3 - bases_left) % 3
                current_phase = Phase(next_phase)
    
    def fix_phase(self, genome_sequence: dict[str, str]) -> bool:
        """
        Fix incorrect phase information in CDS segments.
        
        This method checks if the phase information is consistent with the
        coding sequence and fixes it if necessary.
        
        Args:
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            True if any phase was corrected, False otherwise
        """
        if not self.has_cds:
            return False
        
        # First, infer the phase for any CDS segments with missing phase
        self.infer_phase(genome_sequence)
        
        # Get the coding sequence
        coding_seq = self.cds(genome_sequence)
        
        # Check if the coding sequence length is a multiple of 3
        if len(coding_seq) % 3 != 0:
            # The sequence is not a multiple of 3, so we need to adjust the phase
            # of the last CDS segment to make it a multiple of 3
            if self.strand == Strand.POSITIVE:
                last_cds = max(self.cds_list, key=lambda x: x.end)
            else:
                last_cds = min(self.cds_list, key=lambda x: x.start)
            
            # Calculate how many bases we need to trim
            extra_bases = len(coding_seq) % 3
            
            # Adjust the CDS end position
            if self.strand == Strand.POSITIVE:
                last_cds.end -= extra_bases
            else:
                last_cds.start += extra_bases
            
            # Reset the cached CDS sequence
            self._cds = None
            
            return True
        
        return False
    
    @property
    def introns(self) -> list[Intron]:
        """
        Get the introns for this transcript.
        
        Introns are calculated lazily (only when needed) and cached.
        """
        if self._introns is None:
            self._calc_introns()
        return self._introns if self._introns is not None else []
    
    def _calc_introns(self) -> None:
        """
        Calculate introns from exons.
        
        Introns are the regions between adjacent exons.
        """
        # Initialize empty list
        self._introns = []
        
        # Need at least 2 exons to have introns
        if len(self.exons) < 2:
            return
        
        # Sort exons by position
        sorted_exons = sorted(self.exons, key=lambda x: x.start)
        
        # Create introns between adjacent exons
        for i in range(len(sorted_exons) - 1):
            exon1 = sorted_exons[i]
            exon2 = sorted_exons[i + 1]
            
            # Create intron between these exons
            intron_id = f"{self.id}_intron_{i+1}"
            intron_start = exon1.end
            intron_end = exon2.start
            
            # Only create intron if there's space between exons
            if intron_start < intron_end:
                intron = Intron(
                    id=intron_id,
                    chrom=self.chrom,
                    start=intron_start,
                    end=intron_end,
                    strand=self.strand
                )
                self._introns.append(intron)
    
    def infer_splice_sites(self) -> None:
        """
        Infer splice sites from exon boundaries.
        
        This creates donor (5') and acceptor (3') splice sites at the
        appropriate exon boundaries.
        """
        # Get introns (will calculate if needed)
        introns = self.introns
        
        # Clear existing splice sites
        self.splice_sites = []
        
        # Create splice sites for each intron
        for i, intron in enumerate(introns):
            # Default splice site size (2 bases)
            splice_site_size = 2
            
            # Create donor site (at the start of the intron)
            if self.strand == Strand.POSITIVE:
                # On positive strand, donor is at the end of the preceding exon
                donor_start = intron.start - splice_site_size
                donor_end = intron.start
            else:
                # On negative strand, donor is at the start of the following exon
                donor_start = intron.end
                donor_end = intron.end + splice_site_size
            
            donor_id = f"{self.id}_donor_{i+1}"
            donor_site = SpliceSite(
                id=donor_id,
                chrom=self.chrom,
                start=donor_start,
                end=donor_end,
                strand=self.strand,
                site_type=SpliceSiteType.DONOR
            )
            self.splice_sites.append(donor_site)
            
            # Create acceptor site (at the end of the intron)
            if self.strand == Strand.POSITIVE:
                # On positive strand, acceptor is at the start of the following exon
                acceptor_start = intron.end - splice_site_size
                acceptor_end = intron.end
            else:
                # On negative strand, acceptor is at the end of the preceding exon
                acceptor_start = intron.start
                acceptor_end = intron.start + splice_site_size
            
            acceptor_id = f"{self.id}_acceptor_{i+1}"
            acceptor_site = SpliceSite(
                id=acceptor_id,
                chrom=self.chrom,
                start=acceptor_start,
                end=acceptor_end,
                strand=self.strand,
                site_type=SpliceSiteType.ACCEPTOR
            )
            self.splice_sites.append(acceptor_site)
            
            # Set the donor and acceptor sites on the intron
            intron.donor_site = donor_site
            intron.acceptor_site = acceptor_site
    
    def infer_splice_regions(self, exon_size: int = 3, intron_min: int = 3, intron_max: int = 8) -> None:
        """
        Infer splice regions from exon and intron boundaries.
        
        Splice regions are the areas surrounding splice sites that may contain
        regulatory elements important for splicing.
        
        Args:
            exon_size: Number of bases into the exon to include in the splice region
            intron_min: Minimum number of bases into the intron to include
            intron_max: Maximum number of bases into the intron to include
        """
        # First make sure we have splice sites
        if not self.splice_sites:
            self.infer_splice_sites()
        
        # Create splice regions for each splice site
        for splice_site in self.splice_sites:
            # Create a unique ID for this splice region
            region_id = f"{splice_site.id}_region"
            
            # Determine if this is a donor or acceptor site
            site_type = splice_site.site_type
            
            # Create exonic splice region
            if site_type == SpliceSiteType.DONOR:
                if self.strand == Strand.POSITIVE:
                    # For positive strand, exonic region is before the donor site
                    exon_start = splice_site.start - exon_size
                    exon_end = splice_site.start
                else:
                    # For negative strand, exonic region is after the donor site
                    exon_start = splice_site.end
                    exon_end = splice_site.end + exon_size
            else:  # Acceptor site
                if self.strand == Strand.POSITIVE:
                    # For positive strand, exonic region is after the acceptor site
                    exon_start = splice_site.end
                    exon_end = splice_site.end + exon_size
                else:
                    # For negative strand, exonic region is before the acceptor site
                    exon_start = splice_site.start - exon_size
                    exon_end = splice_site.start
            
            # Create the exonic splice region
            exonic_region = SpliceRegion(
                id=f"{region_id}_exonic",
                chrom=self.chrom,
                start=exon_start,
                end=exon_end,
                strand=self.strand,
                site_type=site_type,
                exonic=True
            )
            
            # Create intronic splice region
            if site_type == SpliceSiteType.DONOR:
                if self.strand == Strand.POSITIVE:
                    # For positive strand, intronic region is after the donor site
                    intron_start = splice_site.end
                    intron_end = splice_site.end + intron_max
                else:
                    # For negative strand, intronic region is before the donor site
                    intron_start = splice_site.start - intron_max
                    intron_end = splice_site.start
            else:  # Acceptor site
                if self.strand == Strand.POSITIVE:
                    # For positive strand, intronic region is before the acceptor site
                    intron_start = splice_site.start - intron_max
                    intron_end = splice_site.start
                else:
                    # For negative strand, intronic region is after the acceptor site
                    intron_start = splice_site.end
                    intron_end = splice_site.end + intron_max
            
            # Create the intronic splice region
            intronic_region = SpliceRegion(
                id=f"{region_id}_intronic",
                chrom=self.chrom,
                start=intron_start,
                end=intron_end,
                strand=self.strand,
                site_type=site_type,
                exonic=False
            )
            
            # Add the splice regions to the transcript
            self.splice_regions.append(exonic_region)
            self.splice_regions.append(intronic_region)
    
    def base_num_to_mrna_pos(self, genomic_pos: int) -> int:
        """
        Calculate the position in mRNA for a genomic position.
        
        Args:
            genomic_pos: Genomic position
            
        Returns:
            Position in mRNA (0-based) or -1 if not in transcript
        """
        # Check if the position is in the transcript
        if not self.intersects_point(genomic_pos):
            return -1
        
        # Sort exons by position based on strand
        if self.strand == Strand.POSITIVE:
            sorted_exons = sorted(self.exons, key=lambda x: x.start)
        else:
            sorted_exons = sorted(self.exons, key=lambda x: x.start, reverse=True)
        
        # Calculate the position in mRNA
        mrna_pos = 0
        for exon in sorted_exons:
            if exon.intersects_point(genomic_pos):
                # Found the exon containing the position
                if self.strand == Strand.POSITIVE:
                    # For positive strand, add the distance from exon start
                    mrna_pos += (genomic_pos - exon.start)
                else:
                    # For negative strand, add the distance from exon end
                    mrna_pos += (exon.end - 1 - genomic_pos)
                return mrna_pos
            else:
                # Add the length of this exon to the mRNA position
                mrna_pos += exon.length
                
                # If we've passed the position, we're done
                if (self.strand == Strand.POSITIVE and exon.start > genomic_pos) or \
                   (self.strand == Strand.NEGATIVE and exon.end <= genomic_pos):
                    break
        
        return -1  # Position not found in any exon
    
    def codon_at(self, genomic_pos: int, genome_sequence: dict[str, str]) -> str:
        """
        Get the codon at a specific genomic position.
        
        Args:
            genomic_pos: Genomic position
            genome_sequence: Dictionary mapping chromosome names to sequences
            
        Returns:
            The codon as a string, or empty string if the position is not in a CDS
        """
        # Convert genomic position to CDS position
        cds_pos = self.genomic_pos_to_cds_pos(genomic_pos, genome_sequence)
        if cds_pos < 0:
            return ""
        
        # Get the coding sequence
        coding_seq = self.cds(genome_sequence)
        
        # Calculate the codon start position
        codon_start = (cds_pos // 3) * 3
        
        # Check if we have enough sequence for a complete codon
        if codon_start + 3 > len(coding_seq):
            return coding_seq[codon_start:]
        
        return coding_seq[codon_start:codon_start + 3]

