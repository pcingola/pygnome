"""
DnaStringArray class for efficiently storing millions of DNA strings in a single array.
"""

import numpy as np
from pygnome.sequences.alphabets import BITS_PER_NUCLEOTIDE, BYTES_PER_INT, DNA_BITS_TO_NT, DNA_NT_TO_BITS, NUCLEOTIDES_PER_INT

DEFAULT_CAPACITY = 1024*1024  # 1 MB default capacity for the data array
DEFAULT_NUMBER_OF_STRINGS = 100 * 1000  # Default number of strings the array can hold
DEFAULT_ENSURE_CAPACITY_STRINGS = 1000  # Default number of sequences to ensure capacity for
DEFAULT_ENSURE_CAPACITY_BYTES = 1000  # Default number of bytes to ensure capacity for

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
        initial_ints = initial_data_bytes // BYTES_PER_INT
        if initial_ints < 1:
            initial_ints = 1
        
        # Initialize data array
        self._data = np.zeros(initial_ints, dtype=np.uint32)
        self._total_nucleotides = 0  # Total number of nucleotides stored

    @property
    def _NT_TO_BITS(self) -> dict[str, int]:
        """Mapping from nucleotide characters to 2-bit values."""
        return DNA_NT_TO_BITS
    
    @property
    def _BITS_TO_NT(self) -> list[str]:
        return DNA_BITS_TO_NT

    def _normalize_nucleotide(self, nt: str) -> str:
        """
        Normalize a nucleotide character.
        """
        return nt.upper()


    def add(self, sequence: str) -> int:
        """
        Add a DNA sequence to the array.
        Returns:
            The index of the added sequence
        """
        # Ensure enough capacity for the new sequence
        self._ensure_capacity(len(sequence), 1)
        
        # Store the starting position of the new sequence
        start_pos = self._positions[self._count]
        
        # Add the sequence to the data array
        for i, nt in enumerate(sequence):
            self._set_nucleotide(start_pos + i, nt)
        # Update the total number of nucleotides
        self._total_nucleotides += len(sequence)
        # Update the end position of the new sequence
        self._positions[self._count + 1] = start_pos + len(sequence)
        # Increment the count of sequences
        self._count += 1
        return self._count - 1  # Return the index of the added sequence
        
    def add_multiple(self, sequences: list[str]) -> list[int]:
        """
        Add multiple DNA sequences to the array.
        
        Args:
            sequences: List of DNA sequences to add
            
        Returns:
            List of indices for the added sequences
        """
        # Calculate total nucleotides to add
        total_nucleotides = sum(len(seq) for seq in sequences)
        
        # Ensure enough capacity for all sequences
        self._ensure_capacity(total_nucleotides, len(sequences))
        
        # Add each sequence and collect indices
        indices = []
        for seq in sequences:
            idx = self.add(seq)
            indices.append(idx)
            
        return indices

    def _ensure_capacity(self, additional_nucleotides: int, additional_sequences: int):
        self._ensure_strings_capacity(additional_sequences)
        self._ensure_data_capacity(additional_nucleotides)


    def _ensure_strings_capacity(self, additional_sequences: int = 1):
        """
        Ensure the positions array has enough capacity for additional sequences and nucleotides.

        Args:
            additional_nucleotides: Number of additional nucleotides to accommodate
            additional_sequences: Number of additional sequences to accommodate
        """
        # Ensure positions array has enough capacity for additional sequences
        required_capacity = self._count + additional_sequences + 1
        if required_capacity > len(self._positions):
            # Double the capacity or grow to required size, whichever is larger
            new_size = max(len(self._positions) * 2, required_capacity)
            new_positions = np.zeros(new_size, dtype=np.uint64)

            # Copy existing positions
            new_positions[:len(self._positions)] = self._positions
            self._positions = new_positions

    
    def _ensure_data_capacity(self, additional_nucleotides: int):
        """
        Ensure the data array has enough capacity for additional nucleotides.
        
        Args:
            additional_nucleotides: Number of additional nucleotides to accommodate
        """
        total_nucleotides = self._total_nucleotides + additional_nucleotides
        required_ints = (total_nucleotides + NUCLEOTIDES_PER_INT - 1) // NUCLEOTIDES_PER_INT
        
        if required_ints > len(self._data):
            # Double the capacity or grow to required size, whichever is larger
            new_size = max(len(self._data) * 2, required_ints)
            new_data = np.zeros(new_size, dtype=np.uint32)
            
            # Copy existing data
            current_ints = (self._total_nucleotides + NUCLEOTIDES_PER_INT - 1) // NUCLEOTIDES_PER_INT
            new_data[:current_ints] = self._data[:current_ints]
            self._data = new_data
    

    def _get_nucleotide(self, nucleotide_idx: int) -> str:
        """
        Get the nucleotide at a specific position.
        This method retrieves the nucleotide at the given index from the packed data.
        """
        int_idx = nucleotide_idx // NUCLEOTIDES_PER_INT
        pos_in_int = nucleotide_idx % NUCLEOTIDES_PER_INT
        value = (self._data[int_idx] >> (pos_in_int * BITS_PER_NUCLEOTIDE)) & 0b11
        return self._BITS_TO_NT[value]


    def get(self, idx: int) -> str:
        """Get a sequence by its index."""
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
        
        # Get start and end positions for this sequence
        start_pos = self._positions[idx]
        end_pos = self._positions[idx + 1]
        length = int(end_pos - start_pos)
        
        # Extract each nucleotide
        dna_str = ''
        for i in range(length):
            dna_str += self._get_nucleotide(start_pos + i)

        return dna_str
        
    def get_subsequence(self, idx: int, start: int, length: int | None = None) -> str:
        """
        Extract a subsequence from a sequence in the array.
        
        Args:
            idx: Index of the sequence
            start: Start position within the sequence (can be negative)
            length: Length of subsequence to extract (default: to end of sequence)
            
        Returns:
            The extracted subsequence
        """
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
            
        seq_start = self._positions[idx]
        seq_end = self._positions[idx + 1]
        seq_length = int(seq_end - seq_start)
        
        # Handle negative start index
        if start < 0:
            start = seq_length + start
            
        if start < 0 or start >= seq_length:
            raise IndexError("Subsequence start index out of range")
            
        if length is None:
            length = seq_length - start
            
        if length < 0 or start + length > seq_length:
            raise IndexError("Subsequence length out of range")
            
        # Extract the subsequence
        result = ''
        for i in range(length):
            result += self._get_nucleotide(seq_start + start + i)
            
        return result

    def get_length(self, idx: int) -> int:
        """Get the length of a sequence."""
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
        return int(self._positions[idx + 1] - self._positions[idx])    

    def __getitem__(self, idx: int) -> str:
        """Get a sequence by its index."""
        return self.get(idx)
    
    def __len__(self) -> int:
        """Return the number of sequences in the array."""
        return self._count
    
    def _set_nucleotide(self, nucleotide_idx: int, value: str):
        """
        Set the nucleotide at a specific index to a new value.
        
        Args:
            nucleotide_idx: Global index of the nucleotide
            value: New nucleotide value (A, C, G, T)
        """
        int_idx = nucleotide_idx // NUCLEOTIDES_PER_INT
        pos_in_int = nucleotide_idx % NUCLEOTIDES_PER_INT
        bit_value = self._NT_TO_BITS.get(value, 0)
        # Clear the bits for this nucleotide
        self._data[int_idx] &= ~(0b11 << (pos_in_int * BITS_PER_NUCLEOTIDE))
        # Set the new bits for this nucleotide
        self._data[int_idx] |= (bit_value << (pos_in_int * BITS_PER_NUCLEOTIDE))


    def trim(self) -> None:
        """
        Trim the internal arrays to their actual used size.
        
        This method reduces memory usage by resizing the data and positions arrays
        to match the actual amount of data stored. This is particularly useful
        before serialization to avoid storing large amounts of unused memory.
        
        Returns:
            None
        """
        if self._count == 0:
            # Reset to minimal arrays if empty
            self._data = np.zeros(1, dtype=np.uint32)
            self._positions = np.zeros(2, dtype=np.uint64)  # Need count+1 positions
            return
            
        # Calculate how many integers we actually need for the data
        required_ints = (self._total_nucleotides + NUCLEOTIDES_PER_INT - 1) // NUCLEOTIDES_PER_INT
        
        # Trim the data array if it's larger than needed
        if len(self._data) > required_ints:
            self._data = self._data[:required_ints].copy()
            
        # Trim the positions array if it's larger than needed
        if len(self._positions) > self._count + 1:
            self._positions = self._positions[:self._count + 1].copy()
        
    def __str__(self) -> str:
        """Return a string representation."""
        count, total_nt, bits_per_nt = self.get_stats()
        return f"DnaStringArray(count={count}, total_nt={total_nt}, bits_per_nt={bits_per_nt:.2f})"
    
    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        return self.__str__()
        
    def get_stats(self) -> tuple[int, int, float]:
        """
        Get statistics about memory usage.
        
        Returns:
            Tuple of (sequence count, total nucleotides, bits per nucleotide)
        """
        if self._count == 0:
            return 0, 0, 0
            
        # Calculate memory usage
        data_bytes = len(self._data) * BYTES_PER_INT
        positions_bytes = len(self._positions) * 8  # 8 bytes per uint64
        total_bytes = data_bytes + positions_bytes
        
        # Calculate bits per nucleotide
        if self._total_nucleotides > 0:
            bits_per_nt = (total_bytes * 8) / self._total_nucleotides
        else:
            bits_per_nt = 0
            
        return self._count, self._total_nucleotides, bits_per_nt
        
    def __getstate__(self):
        """
        Prepare the object for pickling.
        
        This method is called by pickle before serialization.
        It trims the arrays to reduce the serialized size.
        
        Returns:
            The object's state dictionary
        """
        # Trim arrays before pickling
        self.trim()
        
        # Return the object's state
        return self.__dict__
        
    def to_dna_string(self, idx: int) -> 'DnaString':
        """
        Convert a sequence in the array to a DnaString object.
        
        Args:
            idx: Index of the sequence to convert
            
        Returns:
            A DnaString object containing the sequence
        """
        from pygnome.sequences.dna_string import DnaString
        
        if idx < 0 or idx >= self._count:
            raise IndexError("Sequence index out of range")
            
        # Get the sequence as a string
        sequence = self.get(idx)
        
        # Create a new DnaString object
        return DnaString(sequence)