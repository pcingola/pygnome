"""
DnaStringArray class for efficiently storing millions of DNA strings in a single array.
"""

from typing import List, Optional, Tuple
import numpy as np

from .dna_string import DnaString

DEFAULT_CAPACITY = 1024*1024  # 1 MB default capacity for the data array
DEFAULT_NUMBER_OF_STRINGS = 100 * 1000  # Default number of strings the array can hold

class DnaStringArray:
    """
    Efficient storage for millions of small DNA strings in a single NumPy array.
    
    This class stores multiple DNA sequences using the same 2-bit encoding as DnaString
    (A=00, C=01, G=10, T=11), but packs all sequences into a single contiguous array
    for improved memory efficiency when dealing with large numbers of sequences.
    
    The implementation uses:
    - A single uint32 array to store all sequences' data
    - A positions array to track the start position of each sequence
    
    This approach is particularly efficient for storing millions of small sequences
    (like k-mers, reads, or other short DNA fragments) as it:
    1. Minimizes memory overhead compared to individual DnaString objects
    2. Improves cache locality for faster access patterns
    3. Reduces memory fragmentation
    """
    
    # Constants
    NUCLEOTIDES_PER_INT = 16  # Number of nucleotides stored in one integer
    BITS_PER_NUCLEOTIDE = 2   # Number of bits per nucleotide
    BYTES_PER_INT = 4         # Number of bytes in a uint32
    
    # Nucleotide mappings
    _NT_TO_BITS = {
        'A': 0,  # 00
        'C': 1,  # 01
        'G': 2,  # 10
        'T': 3,  # 11
    }
    
    _BITS_TO_NT = {
        0: 'A',  # 00
        1: 'C',  # 01
        2: 'G',  # 10
        3: 'T',  # 11
    }
    
    def __init__(self, initial_data_bytes: int = DEFAULT_CAPACITY, initial_strings: int = DEFAULT_NUMBER_OF_STRINGS):
        """
        Initialize an empty DnaStringArray with the specified capacities.
        
        Args:
            initial_data_bytes: Initial capacity in bytes for the data array
            initial_strings: Initial number of strings the array can hold
        """
        self._count = 0  # Number of sequences currently stored
        self._strings_capacity = initial_strings
        
        # Store only the starting positions (in nucleotides) of each sequence
        # The end of sequence i is the start of sequence i+1
        self._positions = np.zeros(initial_strings + 1, dtype=np.uint64)
        
        # Calculate initial data array size in integers
        initial_ints = initial_data_bytes // self.BYTES_PER_INT
        if initial_ints < 1:
            initial_ints = 1
        
        # Initialize data array
        self._data = np.zeros(initial_ints, dtype=np.uint32)
        self._total_nucleotides = 0  # Total number of nucleotides stored
    
    def _normalize_nucleotide(self, nt: str) -> str:
        """
        Normalize a nucleotide character.
        
        Args:
            nt: A nucleotide character
            
        Returns:
            Normalized nucleotide character
        """
        return nt.upper()
    
    def _ensure_strings_capacity(self, additional_sequences: int = 1):
        """
        Ensure the positions array has enough capacity for additional sequences.
        
        Args:
            additional_sequences: Number of additional sequences to accommodate
        """
        required_capacity = self._count + additional_sequences
        if required_capacity >= self._strings_capacity:  # >= because we store n+1 positions
            # Double the capacity or grow to required size, whichever is larger
            new_capacity = max(self._strings_capacity * 2, required_capacity + 1)
            new_positions = np.zeros(new_capacity + 1, dtype=np.uint64)
            new_positions[:self._count + 1] = self._positions[:self._count + 1]
            self._positions = new_positions
            self._strings_capacity = new_capacity
    
    def _ensure_data_capacity(self, additional_nucleotides: int):
        """
        Ensure the data array has enough capacity for additional nucleotides.
        
        Args:
            additional_nucleotides: Number of additional nucleotides to accommodate
        """
        total_nucleotides = self._total_nucleotides + additional_nucleotides
        required_ints = (total_nucleotides + self.NUCLEOTIDES_PER_INT - 1) // self.NUCLEOTIDES_PER_INT
        
        if required_ints > len(self._data):
            # Double the capacity or grow to required size, whichever is larger
            new_size = max(len(self._data) * 2, required_ints)
            new_data = np.zeros(new_size, dtype=np.uint32)
            
            # Copy existing data
            current_ints = (self._total_nucleotides + self.NUCLEOTIDES_PER_INT - 1) // self.NUCLEOTIDES_PER_INT
            new_data[:current_ints] = self._data[:current_ints]
            self._data = new_data
    
    def add(self, sequence: str) -> int:
        """
        Add a DNA sequence to the array.
        
        Args:
            sequence: A string containing DNA nucleotides (A, C, G, T)
            
        Returns:
            The index of the added sequence
        """
        # Ensure we have capacity for one more sequence
        self._ensure_strings_capacity(1)
        
        length = len(sequence)
        
        # Ensure data array has enough space
        self._ensure_data_capacity(length)
        
        # Store the starting position of this sequence
        seq_idx = self._count
        self._positions[seq_idx] = self._total_nucleotides
        
        # Process sequence in chunks of NUCLEOTIDES_PER_INT nucleotides
        for i in range(0, length, self.NUCLEOTIDES_PER_INT):
            chunk = sequence[i:i+self.NUCLEOTIDES_PER_INT]
            value = 0
            
            # Pack nucleotides into an integer
            for j, nt in enumerate(chunk):
                nt = self._normalize_nucleotide(nt)
                # Use default value (0) for any nucleotide not in our mapping
                bit_value = self._NT_TO_BITS.get(nt, 0)
                # Shift and set the bits for this nucleotide
                value |= (bit_value << (j * self.BITS_PER_NUCLEOTIDE))
            
            # Calculate the correct data index based on the total nucleotides
            data_idx = (self._total_nucleotides + i) // self.NUCLEOTIDES_PER_INT
            self._data[data_idx] = value
        
        # Update total nucleotides and sequence count
        self._total_nucleotides += length
        self._positions[seq_idx + 1] = self._total_nucleotides
        self._count += 1
        
        return seq_idx
    
    def add_multiple(self, sequences: List[str]) -> List[int]:
        """
        Add multiple DNA sequences to the array.
        
        Args:
            sequences: A list of strings containing DNA nucleotides
            
        Returns:
            A list of indices for the added sequences
        """
        # Pre-allocate capacity for all sequences
        self._ensure_strings_capacity(len(sequences))
        
        # Calculate total nucleotides needed
        total_length = sum(len(seq) for seq in sequences)
        
        # Ensure data array has enough space
        self._ensure_data_capacity(total_length)
        
        # Add each sequence
        start_idx = self._count
        for seq in sequences:
            self.add(seq)
        
        # Return indices of added sequences
        return list(range(start_idx, self._count))
    
    def _get_nucleotide_value(self, nucleotide_idx: int) -> int:
        """
        Get the 2-bit value for a nucleotide at a specific position.
        
        Args:
            nucleotide_idx: Global index of the nucleotide
            
        Returns:
            The 2-bit value for the nucleotide
        """
        int_idx = nucleotide_idx // self.NUCLEOTIDES_PER_INT
        pos_in_int = nucleotide_idx % self.NUCLEOTIDES_PER_INT
        return (self._data[int_idx] >> (pos_in_int * self.BITS_PER_NUCLEOTIDE)) & 0b11
    
    def get_sequence(self, idx: int) -> str:
        """
        Get a sequence by its index.
        
        Args:
            idx: Index of the sequence
            
        Returns:
            The sequence as a string
        """
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
        
        # Get start and end positions for this sequence
        start_pos = self._positions[idx]
        end_pos = self._positions[idx + 1]
        length = int(end_pos - start_pos)
        
        # Special case for empty sequence
        if length == 0:
            return ""
        
        result = []
        
        # Extract each nucleotide
        for i in range(length):
            nucleotide_idx = start_pos + i
            int_idx = nucleotide_idx // self.NUCLEOTIDES_PER_INT
            pos_in_int = nucleotide_idx % self.NUCLEOTIDES_PER_INT
            
            # Extract the bits for this nucleotide
            value = (self._data[int_idx] >> (pos_in_int * self.BITS_PER_NUCLEOTIDE)) & 0b11
            result.append(self._BITS_TO_NT[value])
        
        return ''.join(result)
    
    def __getitem__(self, idx: int) -> str:
        """
        Get a sequence by its index.
        
        Args:
            idx: Index of the sequence
            
        Returns:
            The sequence as a string
        """
        return self.get_sequence(idx)
    
    def __len__(self) -> int:
        """Return the number of sequences in the array."""
        return self._count
    
    def get_length(self, idx: int) -> int:
        """
        Get the length of a sequence.
        
        Args:
            idx: Index of the sequence
            
        Returns:
            Length of the sequence in nucleotides
        """
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
        
        return int(self._positions[idx + 1] - self._positions[idx])
    
    def get_subsequence(self, idx: int, start: int, length: Optional[int] = None) -> str:
        """
        Get a subsequence from a sequence.
        
        Args:
            idx: Index of the sequence
            start: Starting position within the sequence (0-based)
            length: Length of subsequence to extract, or None for the rest of the sequence
            
        Returns:
            The extracted subsequence
        """
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
        
        # Get start and end positions for this sequence
        seq_start = self._positions[idx]
        seq_end = self._positions[idx + 1]
        seq_length = int(seq_end - seq_start)  # Convert to int to handle negative indices
        
        # Handle negative start index
        if start < 0:
            start = seq_length + start  # Convert to positive index
        
        if start < 0 or start >= seq_length:
            raise IndexError("Subsequence start index out of range")
        
        if length is None:
            length = seq_length - start
        
        if length < 0 or start + length > seq_length:
            raise IndexError("Subsequence length out of range")
        
        result = []
        
        # Extract each nucleotide in the subsequence
        for i in range(length):
            nucleotide_idx = seq_start + start + i
            value = self._get_nucleotide_value(nucleotide_idx)
            result.append(self._BITS_TO_NT[value])
        
        return ''.join(result)
    
    def to_dna_string(self, idx: int) -> DnaString:
        """
        Convert a sequence to a DnaString object.
        
        Args:
            idx: Index of the sequence
            
        Returns:
            A DnaString object representing the sequence
        """
        return DnaString(self.get_sequence(idx))
    
    def get_stats(self) -> Tuple[int, int, float]:
        """
        Get statistics about memory usage.
        
        Returns:
            A tuple containing:
            - Number of sequences
            - Total number of nucleotides
            - Average bits per nucleotide
        """
        total_nucleotides = self._total_nucleotides
        data_ints_used = (total_nucleotides + self.NUCLEOTIDES_PER_INT - 1) // self.NUCLEOTIDES_PER_INT
        data_bits = data_ints_used * 32  # 32 bits per uint32
        
        # Add bits used for positions array
        positions_bits = (self._count + 1) * 64  # 64 bits per uint64
        total_bits = data_bits + positions_bits
        
        if total_nucleotides == 0:
            avg_bits_per_nt = 0
        else:
            avg_bits_per_nt = total_bits / total_nucleotides
        
        return self._count, total_nucleotides, avg_bits_per_nt
    
    def __str__(self) -> str:
        """Return a string representation."""
        count, total_nt, bits_per_nt = self.get_stats()
        return f"DnaStringArray(count={count}, total_nt={total_nt}, bits_per_nt={bits_per_nt:.2f})"
    
    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        return self.__str__()