"""
Base class for efficient 2-bit representation of nucleotide sequences.
"""

from abc import ABC, abstractmethod
import numpy as np

from pygnome.sequences.alphabets import BITS_PER_NUCLEOTIDE, NUCLEOTIDES_PER_INT, PREVIEW_LENGTH


class BaseSequence(ABC):
    """
    Base class for efficient 2-bit representation of nucleotide sequences.
    
    This class provides common functionality for storing nucleotide sequences
    using 2 bits per nucleotide, allowing 16 nucleotides to be packed into a
    single 32-bit integer.
    """
    
    @property
    @abstractmethod
    def _NT_TO_BITS(self) -> dict[str, int]:
        """Mapping from nucleotide characters to 2-bit values."""
        pass
    
    @property
    @abstractmethod
    def _BITS_TO_NT(self) -> dict[int, str]:
        """Mapping from 2-bit values to nucleotide characters."""
        pass
    
    def __init__(self, sequence: str):
        """
        Initialize a sequence from a string.
        
        Args:
            sequence: A string containing nucleotides.
        """
        self.length = len(sequence)
        
        # Calculate how many integers we need to store the sequence
        num_ints = (self.length + NUCLEOTIDES_PER_INT - 1) // NUCLEOTIDES_PER_INT
        self._data = np.zeros(num_ints, dtype=np.uint32)
        
        # Process sequence in chunks of NUCLEOTIDES_PER_INT nucleotides
        for i in range(0, self.length, NUCLEOTIDES_PER_INT):
            chunk = sequence[i:i+NUCLEOTIDES_PER_INT]
            value = 0
            
            # Pack nucleotides into an integer
            for j, nt in enumerate(chunk):
                nt = self._normalize_nucleotide(nt)
                # Use default value (0) for any nucleotide not in our mapping
                bit_value = self._NT_TO_BITS.get(nt, 0)
                # Shift and set the bits for this nucleotide
                value |= (bit_value << (j * BITS_PER_NUCLEOTIDE))
            
            self._data[i // NUCLEOTIDES_PER_INT] = value
    
    def _normalize_nucleotide(self, nt: str) -> str:
        """
        Normalize a nucleotide character.
        """
        return nt.upper()
    
    def __len__(self) -> int:
        """Return the length of the sequence."""
        return self.length
    
    def __str__(self) -> str:
        """Return the sequence as a string."""
        return self.to_string()
    
    def _get_nucleotide_value(self, index: int) -> int:
        """
        Get the 2-bit value for the nucleotide at the given index.
        """
        # Calculate which integer contains this nucleotide
        int_idx = index // NUCLEOTIDES_PER_INT
        # Calculate position within the integer
        pos_in_int = index % NUCLEOTIDES_PER_INT
        # Extract the bits for this nucleotide
        return (self._data[int_idx] >> (pos_in_int * BITS_PER_NUCLEOTIDE)) & 0b11
    
    def __getitem__(self, key) -> str:
        """
        Get a nucleotide or subsequence.
        """
        if isinstance(key, int):
            # Handle negative indices
            if key < 0:
                key += self.length
            
            if key < 0 or key >= self.length:
                raise IndexError(f"{self._class_name} index out of range")
            
            # Get the 2-bit value for this nucleotide
            value = self._get_nucleotide_value(key)
            
            return self._BITS_TO_NT[value]
        
        elif isinstance(key, slice):
            # Handle slicing
            start, stop, step = key.indices(self.length)
            
            if step == 1:
                # Optimize for continuous slices
                return self.substring(start, stop - start)
            else:
                # Handle step != 1
                result = ""
                for i in range(start, stop, step):
                    result += self[i]
                return result
        
        else:
            raise TypeError(f"{self._class_name} indices must be integers or slices")
    
    def to_string(self) -> str:
        return self.substring(0, self.length)
    
    def substring(self, start: int, length: int |  None = None) -> str:
        """
        Extract a substring from the sequence.
        """
        if start < 0:
            start += self.length
        
        if start < 0 or start >= self.length:
            raise IndexError("Substring start index out of range")
        
        if length is None:
            length = self.length - start
        
        if length < 0 or start + length > self.length:
            raise IndexError("Substring length out of range")
        
        result = ""
        
        for i in range(start, start + length):
            # Get the 2-bit value for this nucleotide
            value = self._get_nucleotide_value(i)
            
            result += self._BITS_TO_NT[value]
        
        return result
    
    def __eq__(self, other) -> bool:
        """Check if two sequence objects are equal."""
        if not isinstance(other, self.__class__):
            return False
        
        if self.length != other.length:
            return False
        
        return np.array_equal(self._data, other._data)
    
    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        return self.to_string()
