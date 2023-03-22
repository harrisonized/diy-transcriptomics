# FastQC

# cd diy-transcriptomics
conda activate rnaseq
fastqc data/fastq/SRR8668755_1M_subsample.fastq -t 8
# outputs SRR8668755_1M_subsample_fastqc.html in the same directory


# Kallisto

# indexing
kallisto index --index='data/lecture_2/ref/fasta/Homo_sapiens.GRCh38.cdna.all.index' data/lecture_2/ref/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz

# [build]loading fasta file data/ref/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz
# [build]k-mer length: 31
# [build]warning: clipped off poly-A tail (longer than 10)
#       from 1517 target sequences
# [build]warning: replaced 100005 non-ACGUT characters in the input sequence
#       with pseudorandom nucleotides
# [build]counting k-mers ... done.
# [build]building target de Bruijn graph ...  done 
# [build]creating equivalence classes ...  done
# [build]target de Bruijn graph has 1229410 contigs and contains 116457368 k-mers 


# quantification (mapping)
# index, output, threads, single, single vs. paired end, avg length and stdev for library prep
kallisto quant --index='data/lecture_2/ref/fasta/Homo_sapiens.GRCh38.cdna.all.index' --output-dir=data/lecture_2/mapped_reads/SRR8668755_subset -t 8 --single -l 250 -s 30 data/lecture_2/fastq/SRR8668755_1M_subsample.fastq.gz &> logs/SRR8668755_subset.log
# outputs abundance.h5, abundance.tsv, run_info.json


# multiqc
# generate multiqc_report.html (summary report)
multiqc -d data
#  /// MultiQC 🔍 | v1.14
#
# |           multiqc | Prepending directory to sample names
# |           multiqc | Search path : /home/harrisonized/github/diy-transcriptomics/data
# |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 50/50  
# |            fastqc | Found 1 reports
# |           multiqc | Compressing plot data
# |           multiqc | Report      : multiqc_report.html
# |           multiqc | Data        : multiqc_data
# |           multiqc | MultiQC complete

