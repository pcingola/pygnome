"""Utility functions for DNA and RNA sequence operations."""

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    # Define the complement mapping
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    
    # Reverse the sequence and get the complement of each nucleotide
    rev_comp = ''.join(complement.get(nt, nt) for nt in reversed(sequence))
    
    return rev_comp

def dna_to_rna(sequence: str) -> str:
    """Convert a DNA sequence to RNA by replacing T with U."""
    return sequence.replace('T', 'U').replace('t', 'u')