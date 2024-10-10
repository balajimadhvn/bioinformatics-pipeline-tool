import subprocess
import os
import gzip
from Bio import SeqIO
import argparse

print("""
$$$$$$\\                      $$$$$$\\                                               $$\\     $$$$$$$\\  $$\\   $$\\  $$$$$$\\  
$$  __$$\\                     \\_$$  _|                                              $$ |    $$  __$$\\ $$$\\  $$ |$$  __$$\\ 
$$ /  \\__| $$$$$$\\   $$$$$$\\    $$ |  $$$$$$\\$$$$\\   $$$$$$\\   $$$$$$\\   $$$$$$$\\ $$$$$$\\   $$ |  $$ |$$$$\\ $$ |$$ /  $$ |
\\$$$$$$\\  $$  __$$\\ $$  __$$\\   $$ |  $$  _$$  _$$\\ $$  __$$\\  \\____$$\\ $$  _____|\\_$$  _|  $$$$$$$  |$$ $$\\$$ |$$$$$$$$ |
 \\____$$\\ $$$$$$$$ |$$ /  $$ |  $$ |  $$ / $$ / $$ |$$ /  $$ | $$$$$$$ |$$ /        $$ |    $$  __$$< $$ \\$$$$ |$$  __$$ |
$$\\   $$ |$$   ____|$$ |  $$ |  $$ |  $$ | $$ | $$ |$$ |  $$ |$$  __$$ |$$ |        $$ |$$\\ $$ |  $$ |$$ |\\$$$ |$$ |  $$ |
\\$$$$$$  |\\$$$$$$$\\ \\$$$$$$$ |$$$$$$\\ $$ | $$ | $$ |$$$$$$$  |\\$$$$$$$ |\\$$$$$$$\\   \\$$$$  |$$ |  $$ |$$ | \\$$ |$$ |  $$ |
 \\______/  \\_______| \\____$$ |\\______|\\__| \\__| \\__|$$  ____/  \\_______| \\_______|   \\____/ \\__|  \\__|\\__|  \\__|\\__|  \\__|
                          $$ |                      $$ |                                                                  
                          $$ |                      $$ |                                                                  
                          \\__|                      \\__|                                                                  """)

def run_fastqc(input_file, output_dir):
    """
    Runs FastQC on input FASTQ file.
    
    Parameters:
    input_file (str): Path to input FASTQ file.
    output_dir (str): Output directory for FastQC results.
    """
    fastqc_cmd = ["fastqc", "--outdir", output_dir, input_file]
    subprocess.run(fastqc_cmd, check=True)

def trim_reads(input_file, output_file, min_quality):
    """
    Trims reads from input FASTQ file based on minimum average quality score.
    
    Parameters:
    input_file (str): Path to input FASTQ file.
    output_file (str): Path to output trimmed FASTQ file.
    min_quality (int): Minimum average quality score threshold.
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            # Calculate average quality score for the record
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)

            # Trim based on average quality score
            if avg_quality >= min_quality:
                SeqIO.write(record, out_handle, "fastq")

def align_to_reference(forward_fastq, reverse_fastq, reference_index_prefix, output_bam):
    """
    Aligns paired-end reads to a reference genome using bowtie2 and converts to BAM.
    
    Parameters:
    forward_fastq (str): Path to forward FASTQ file.
    reverse_fastq (str): Path to reverse FASTQ file.
    reference_index_prefix (str): Prefix of the bowtie2 index files for the reference genome.
    output_bam (str): Path to output BAM file after alignment.
    """
    bowtie2_cmd = [
        "bowtie2",
        "-x", reference_index_prefix,
        "-1", forward_fastq,
        "-2", reverse_fastq,
        "|",
        "samtools", "view", "-bS", "-",
        "|",
        "samtools", "sort", "-o", output_bam, "-"
    ]
    subprocess.run(" ".join(bowtie2_cmd), shell=True, check=True)

def remove_duplicates(input_bam, output_bam):
    """
    Removes duplicates from BAM file using samtools.
    
    Parameters:
    input_bam (str): Path to input BAM file.
    output_bam (str): Path to output BAM file after removing duplicates.
    """
    samtools_rmdup_cmd = [
        "samtools", "rmdup", "-s", input_bam, output_bam
    ]
    subprocess.run(samtools_rmdup_cmd, check=True)

def base_quality_recalibration(input_bam, reference_fasta, output_recal_table, output_recalibrated_bam):
    """
    Performs base quality score recalibration (BQSR) using GATK.
    
    Parameters:
    input_bam (str): Path to input BAM file.
    reference_fasta (str): Path to reference genome FASTA file.
    output_recal_table (str): Path to output recalibration table.
    output_recalibrated_bam (str): Path to output recalibrated BAM file.
    """
    gatk_cmd = [
        "gatk", "BaseRecalibrator",
        "-I", input_bam,
        "-R", reference_fasta,
        "--known-sites", "known_sites.vcf",  # Replace with known variant sites (e.g., dbSNP)
        "-O", output_recal_table
    ]
    subprocess.run(gatk_cmd, check=True)
    
    gatk_apply_cmd = [
        "gatk", "ApplyBQSR",
        "-I", input_bam,
        "-R", reference_fasta,
        "--bqsr-recal-file", output_recal_table,
        "-O", output_recalibrated_bam
    ]
    subprocess.run(gatk_apply_cmd, check=True)

def index_bam(input_bam):
    """
    Indexes the BAM file using samtools index.
    
    Parameters:
    input_bam (str): Path to input BAM file.
    """
    samtools_index_cmd = [
        "samtools", "index", input_bam
    ]
    subprocess.run(samtools_index_cmd, check=True)

def count_features(input_bam, annotation_file, output_counts):
    """
    Counts reads mapped to features (e.g., genes) using featureCounts.
    
    Parameters:
    input_bam (str): Path to input BAM file.
    annotation_file (str): Path to annotation GTF/GFF file.
    output_counts (str): Path to output counts file.
    """
    featureCounts_cmd = [
        "featureCounts",
        "-T", "4",  # Number of threads
        "-Q", "30",  # Minimum mapping quality
        "-a", annotation_file,
        "-o", output_counts,
        input_bam
    ]
    subprocess.run(featureCounts_cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess and analyze paired-end sequencing data with quality checks, trimming, alignment, and feature counting."
    )
    
    parser.add_argument(
        "--input_fastq_forward",
        required=True,
        help="Path to input forward FASTQ file."
    )
    parser.add_argument(
        "--input_fastq_reverse",
        required=True,
        help="Path to input reverse FASTQ file."
    )
    parser.add_argument(
        "--min_quality",
        type=int,
        default=20,
        help="Minimum average quality score threshold for trimming reads. Default is 20."
    )
    parser.add_argument(
        "--output_dir",
        default=".",
        help="Directory where FastQC results will be saved and all output files will be placed. Default is the current directory."
    )
    parser.add_argument(
        "--reference_index_prefix",
        required=True,
        help="Prefix of the bowtie2 index files for the reference genome."
    )
    parser.add_argument(
        "--reference_fasta",
        required=True,
        help="Path to the reference genome FASTA file."
    )
    parser.add_argument(
        "--annotation_file",
        required=True,
        help="Path to the annotation GTF/GFF file for featureCounts."
    )
    
    args = parser.parse_args()

    # Ensure output directory exists or create it if not
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Check quality of input FASTQ files using FastQC
    print("Step 1: Running FastQC on input FASTQ files...")
    run_fastqc(args.input_fastq_forward, args.output_dir)
    run_fastqc(args.input_fastq_reverse, args.output_dir)

    # Step 2: Trim reads based on quality for forward and reverse reads
    print("Step 2: Trimming reads based on quality...")
    trimmed_fastq_forward = os.path.join(args.output_dir, "trimmed_forward.fastq")
    trimmed_fastq_reverse = os.path.join(args.output_dir, "trimmed_reverse.fastq")
    trim_reads(args.input_fastq_forward, trimmed_fastq_forward, args.min_quality)
    trim_reads(args.input_fastq_reverse, trimmed_fastq_reverse, args.min_quality)
    
    # Step 3: Check quality of trimmed FASTQ files using FastQC
    print("Step 3: Running FastQC on trimmed FASTQ files...")
    run_fastqc(trimmed_fastq_forward, args.output_dir)
    run_fastqc(trimmed_fastq_reverse, args.output_dir)

    # Step 4: Align trimmed reads to the reference genome and convert to BAM
    print("Step 4: Aligning trimmed reads to the reference genome and converting to BAM...")
    aligned_bam = os.path.join(args.output_dir, "aligned_reads.bam")
    align_to_reference(trimmed_fastq_forward, trimmed_fastq_reverse, args.reference_index_prefix, aligned_bam)

    # Step 5: Remove duplicates from aligned BAM file
    print("Step 5: Removing duplicates from aligned BAM file...")
    dedup_bam = os.path.join(args.output_dir, "deduplicated_reads.bam")
    remove_duplicates(aligned_bam, dedup_bam)

    # Step 6: Perform base quality score recalibration (BQSR)
    print("Step 6: Performing base quality score recalibration (BQSR)...")
    recal_table = os.path.join(args.output_dir, "recalibration_table.txt")
    recalibrated_bam = os.path.join(args.output_dir, "recalibrated_reads.bam")
    base_quality_recalibration(dedup_bam, args.reference_fasta, recal_table, recalibrated_bam)

    # Step 7: Index the recalibrated BAM file
    print("Step 7: Indexing the recalibrated BAM file...")
    index_bam(recalibrated_bam)

    # Step 8: Count reads mapped to features using featureCounts
    print("Step 8: Counting reads mapped to features using featureCounts...")
    feature_counts = os.path.join(args.output_dir, "feature_counts.txt")
    count_features(recalibrated_bam, args.annotation_file, feature_counts)

    print("Processing complete.")

