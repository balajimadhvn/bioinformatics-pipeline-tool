# SeqImpactRNA

SeqImpactRNA is a Python tool for analyzing RNA-seq data, focusing on identifying differentially expressed genes (DEGs) and exploring their impact.

## Usage

- Ensure you have the necessary dependencies installed (`Bio`, `samtools`, `bowtie2`, `gatk`, `featureCounts`).
- Modify the script arguments as per your data and reference genome.
- Run the tool using Python:

  ```bash
  python SeqImpactRNA.py --input_fastq_forward data/forward.fastq --input_fastq_reverse data/reverse.fastq --reference_index_prefix ref_genome/index --reference_fasta ref_genome/genome.fa --annotation_file annotations/genes.gtf
