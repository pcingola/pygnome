"""
Base class for efficient 2-bit representation of nucleotide sequences.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Optional


class BaseSequence(ABC):
    """
    Base class for efficient 2-bit representation of nucleotide sequences.
    
    This class provides common functionality for storing nucleotide sequences
    using 2 bits per nucleotide, allowing 16 nucleotides to be packed into a
    single 32-bit integer.
    """
    
    # Constants
    NUCLEOTIDES_PER_INT = 16  # Number of nucleotides stored in one integer
    BITS_PER_NUCLEOTIDE = 2   # Number of bits per nucleotide
    PREVIEW_LENGTH = 30       # Length of sequence preview in __repr__
    
    @property
    @abstractmethod
    def _NT_TO_BITS(self) -> Dict[str, int]:
        """Mapping from nucleotide characters to 2-bit values."""
        pass
    
    @property
    @abstractmethod
    def _BITS_TO_NT(self) -> Dict[int, str]:
        """Mapping from 2-bit values to nucleotide characters."""
        pass
    
    @property
    @abstractmethod
    def _class_name(self) -> str:
        """Name of the class for use in __repr__."""
        pass
    
    def __init__(self, sequence: str):
        """
        Initialize a sequence from a string.
        
        Args:
            sequence: A string containing nucleotides.
        """
        self.length = len(sequence)
        self._data: List[int] = []
        
        # Process sequence in chunks of NUCLEOTIDES_PER_INT nucleotides
        for i in range(0, self.length, self.NUCLEOTIDES_PER_INT):
            chunk = sequence[i:i+self.NUCLEOTIDES_PER_INT]
            value = 0
            
            # Pack nucleotides into an integer
            for j, nt in enumerate(chunk):
                nt = self._normalize_nucleotide(nt)
                # Use default value (0) for any nucleotide not in our mapping
                bit_value = self._NT_TO_BITS.get(nt, 0)
                # Shift and set the bits for this nucleotide
                value |= (bit_value << (j * self.BITS_PER_NUCLEOTIDE))
            
            self._data.append(value)
    
    def _normalize_nucleotide(self, nt: str) -> str:
        """
        Normalize a nucleotide character.
        
        Args:
            nt: A nucleotide character
            
        Returns:
            Normalized nucleotide character
        """
        return nt.upper()
    
    def __len__(self) -> int:
        """Return the length of the sequence."""
        return self.length
    
    def __str__(self) -> str:
        """Return the sequence as a string."""
        return self.to_string()
    
    def __getitem__(self, key) -> str:
        """
        Get a nucleotide or subsequence.
        
        Args:
            key: An index or slice
            
        Returns:
            A single nucleotide character or a substring
        """
        if isinstance(key, int):
            # Handle negative indices
            if key < 0:
                key += self.length
            
            if key < 0 or key >= self.length:
                raise IndexError(f"{self._class_name} index out of range")
            
            # Calculate which integer contains this nucleotide
            int_idx = key // self.NUCLEOTIDES_PER_INT
            # Calculate position within the integer
            pos_in_int = key % self.NUCLEOTIDES_PER_INT
            # Extract the bits for this nucleotide
            value = (self._data[int_idx] >> (pos_in_int * self.BITS_PER_NUCLEOTIDE)) & 0b11
            
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
        """
        Convert the entire sequence to a string.
        
        Returns:
            The complete sequence as a string
        """
        return self.substring(0, self.length)
    
    def substring(self, start: int, length: Optional[int] = None) -> str:
        """
        Extract a substring from the sequence.
        
        Args:
            start: Starting position (0-based)
            length: Length of substring to extract, or None for the rest of the sequence
            
        Returns:
            The extracted substring
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
            # Calculate which integer contains this nucleotide
            int_idx = i // self.NUCLEOTIDES_PER_INT
            # Calculate position within the integer
            pos_in_int = i % self.NUCLEOTIDES_PER_INT
            # Extract the bits for this nucleotide
            value = (self._data[int_idx] >> (pos_in_int * self.BITS_PER_NUCLEOTIDE)) & 0b11
            
            result += self._BITS_TO_NT[value]
        
        return result
    
    def __eq__(self, other) -> bool:
        """Check if two sequence objects are equal."""
        if not isinstance(other, self.__class__):
            return False
        
        if self.length != other.length:
            return False
        
        return self._data == other._data
    
    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        if self.length <= 2 * self.PREVIEW_LENGTH:
            return f"{self._class_name}('{self.to_string()}')"
        else:
            return f"{self._class_name}('{self.substring(0, self.PREVIEW_LENGTH)}...{self.substring(self.length-self.PREVIEW_LENGTH, self.PREVIEW_LENGTH)}')"