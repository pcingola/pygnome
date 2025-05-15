# Constants
NUCLEOTIDES_PER_INT = 16  # Number of nucleotides stored in one integer
BITS_PER_NUCLEOTIDE = 2   # Number of bits per nucleotide
PREVIEW_LENGTH = 30       # Length of sequence preview in __repr__
BYTES_PER_INT = 4  # Number of bytes in an integer

# Nucleotide mappings
DNA_NT_TO_BITS = {
    'A': 0,  # 00
    'C': 1,  # 01
    'G': 2,  # 10
    'T': 3,  # 11
}

DNA_BITS_TO_NT = [
    'A',  # 00
    'C',  # 01
    'G',  # 10
    'T',  # 11
]


RNA_NT_TO_BITS = {
    'A': 0,  # 00
    'C': 1,  # 01
    'G': 2,  # 10
    'U': 3,  # 11
}

RNA_BITS_TO_NT = [
    'A',  # 00
    'C',  # 01
    'G',  # 10
    'U',  # 11
]