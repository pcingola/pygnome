"""CodonTable class for translating codons to amino acids."""

from enum import Enum
from typing import Dict, Set


class CodonTableType(str, Enum):
    """Enumeration of standard codon table types."""
    STANDARD = "Standard"
    VERTEBRATE_MITOCHONDRIAL = "Vertebrate_Mitochondrial"
    YEAST_MITOCHONDRIAL = "Yeast_Mitochondrial"
    MOLD_MITOCHONDRIAL = "Mold_Mitochondrial"
    PROTOZOAN_MITOCHONDRIAL = "Protozoan_Mitochondrial"
    COELENTERATE = "Coelenterate"
    MITOCHONDRIAL = "Mitochondrial"
    MYCOPLASMA = "Mycoplasma"
    SPIROPLASMA = "Spiroplasma"
    INVERTEBRATE_MITOCHONDRIAL = "Invertebrate_Mitochondrial"
    CILIATE_NUCLEAR = "Ciliate_Nuclear"
    DASYCLADACEAN_NUCLEAR = "Dasycladacean_Nuclear"
    HEXAMITA_NUCLEAR = "Hexamita_Nuclear"
    ECHINODERM_MITOCHONDRIAL = "Echinoderm_Mitochondrial"
    FLATWORM_MITOCHONDRIAL = "Flatworm_Mitochondrial"
    EUPLOTID_NUCLEAR = "Euplotid_Nuclear"
    BACTERIAL_AND_PLANT_PLASTID = "Bacterial_and_Plant_Plastid"
    ALTERNATIVE_YEAST_NUCLEAR = "Alternative_Yeast_Nuclear"
    ASCIDIAN_MITOCHONDRIAL = "Ascidian_Mitochondrial"
    ALTERNATIVE_FLATWORM_MITOCHONDRIAL = "Alternative_Flatworm_Mitochondrial"
    BLEPHARISMA_MACRONUCLEAR = "Blepharisma_Macronuclear"
    CHLOROPHYCEAN_MITOCHONDRIAL = "Chlorophycean_Mitochondrial"
    TREMATODE_MITOCHONDRIAL = "Trematode_Mitochondrial"
    SCENEDESMUS_OBLIQUUS_MITOCHONDRIAL = "Scenedesmus_obliquus_Mitochondrial"
    THRAUSTOCHYTRIUM_MITOCHONDRIAL = "Thraustochytrium_Mitochondrial"


class CodonTable:
    """Translates codons to amino acids using different codon tables."""
    
    # Define standard codon tables
    _TABLES = {
        CodonTableType.STANDARD: {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        },
        CodonTableType.VERTEBRATE_MITOCHONDRIAL: {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "*", "AGG": "*",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        },
        # Additional tables would be defined here
    }
    
    # Define start codons for each table
    _START_CODONS = {
        CodonTableType.STANDARD: {"TTG", "CTG", "ATG"},
        CodonTableType.VERTEBRATE_MITOCHONDRIAL: {"ATT", "ATC", "ATA", "ATG", "GTG"},
        # Additional start codons would be defined here
    }
    
    def __init__(self, table_type: CodonTableType = CodonTableType.STANDARD):
        """Initialize a codon table."""
        self.table_type = table_type
        
        # Initialize the table and start codons
        if table_type in self._TABLES:
            self.table = self._TABLES[table_type]
        else:
            raise ValueError(f"Unsupported codon table type: {table_type}")
        
        if table_type in self._START_CODONS:
            self.start_codons = self._START_CODONS[table_type]
        else:
            self.start_codons = set()
    
    def translate_codon(self, codon: str) -> str:
        """Translate a single codon to an amino acid."""
        codon = codon.upper()
        if len(codon) != 3:
            raise ValueError(f"Invalid codon length: {len(codon)}, must be 3")
        
        return self.table.get(codon, "X")  # Return 'X' for unknown codons
    
    def translate_sequence(self, sequence: str) -> str:
        """Translate a DNA sequence to an amino acid sequence."""
        if len(sequence) % 3 != 0:
            raise ValueError(f"Sequence length ({len(sequence)}) is not a multiple of 3")
        
        amino_acids = []
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            amino_acids.append(self.translate_codon(codon))
        
        return "".join(amino_acids)
    
    def is_start_codon(self, codon: str) -> bool:
        """Check if a codon is a start codon in this table."""
        return codon.upper() in self.start_codons
    
    def is_stop_codon(self, codon: str) -> bool:
        """Check if a codon is a stop codon in this table."""
        return self.translate_codon(codon) == "*"
    
    @classmethod
    def parse_table_definition(cls, definition: str) -> tuple[Dict[str, str], Set[str]]:
        """Parse a codon table definition string."""
        table = {}
        start_codons = set()
        
        # Split the definition into codon/amino acid pairs
        pairs = [pair.strip() for pair in definition.split(",")]
        
        for pair in pairs:
            # Parse each codon/amino acid pair
            codon_aa = pair.strip().split("/")
            if len(codon_aa) != 2:
                continue
                
            codon, aa = codon_aa
            codon = codon.strip()
            aa = aa.strip()
            
            # Check for start codons (marked with +)
            if aa.endswith("+"):
                aa = aa[:-1]
                start_codons.add(codon)
                
            table[codon] = aa
            
        return table, start_codons
    
    @classmethod
    def from_definition(cls, definition: str):
        """Create a codon table from a definition string."""
        table, start_codons = cls.parse_table_definition(definition)
        
        # Create a new instance with a custom table
        instance = cls(CodonTableType.STANDARD)
        instance.table = table
        instance.start_codons = start_codons
        
        return instance
    
    def __str__(self) -> str:
        """Return a string representation of the codon table."""
        return f"CodonTable({self.table_type})"