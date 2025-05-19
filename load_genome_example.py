#!/usr/bin/env python3
"""
Script to load genomic structure and chromosome sequence from GTF and FASTA files.
"""

from pathlib import Path

from pygnome.parsers.genome_loader import GenomeLoader


def main():
    """Load genomic structure and chromosome sequence."""
    # File paths
    gtf_file = Path("data/genomes/GRCh38.mane.1.2.ensembl.chr21.gtf")
    fasta_file = Path("data/genomes/chr21.fa.gz")
    
    # Create a genome loader with verbose mode enabled
    loader = GenomeLoader(
        annotation_file=gtf_file,
        sequence_file=fasta_file,
        genome_name="GRCh38.mane.1.2",
        species="Homo sapiens",
        verbose=True
    )
    
    # Load the genome
    genome = loader.load()
    
    # Print information about the loaded genome
    print(f"Loaded genome: {genome}")
    
    # Print information about chromosomes
    for chrom_name, chrom in genome.chromosomes.items():
        print(f"Chromosome: {chrom}")
        print(f"  Genes: {len(chrom.genes)}")
        
        # Count genes by biotype
        biotype_counts = {}
        for gene in chrom.genes.values():
            biotype = gene.biotype.value if gene.biotype else "unknown"
            biotype_counts[biotype] = biotype_counts.get(biotype, 0) + 1
        
        # Print biotype counts
        for biotype, count in sorted(biotype_counts.items()):
            print(f"    {biotype}: {count} genes")
        
        # Count transcripts
        transcript_count = sum(len(gene.transcripts) for gene in chrom.genes.values())
        print(f"  Transcripts: {transcript_count}")
        
        # Count exons and CDS
        exon_count = 0
        cds_count = 0
        for gene in chrom.genes.values():
            for transcript in gene.transcripts:
                exon_count += len(transcript.exons)
                cds_count += len(transcript.cds_list)
        
        print(f"  Exons: {exon_count}")
        print(f"  CDS: {cds_count}")
    
    # Print a few example genes
    print("\nExample genes:")
    for i, gene in enumerate(sorted(genome.genes.values(), key=lambda g: g.start)[:5]):
        print(f"  {gene.id} ({gene.name}): {gene.chrom}:{gene.start}-{gene.end} ({gene.strand})")
        print(f"    Biotype: {gene.biotype.value if gene.biotype else 'unknown'}")
        print(f"    Transcripts: {len(gene.transcripts)}")
        
        # Print a couple of transcripts for each gene
        for j, transcript in enumerate(gene.transcripts[:2]):
            print(f"      {transcript.id}: {transcript.chrom}:{transcript.start}-{transcript.end}")
            print(f"        Exons: {len(transcript.exons)}, CDS: {len(transcript.cds_list)}")
            
            # If there are more transcripts, indicate how many more
            if j == 1 and len(gene.transcripts) > 2:
                print(f"        ... and {len(gene.transcripts) - 2} more transcripts")


if __name__ == "__main__":
    main()