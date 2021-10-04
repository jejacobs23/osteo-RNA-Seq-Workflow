# osteo-RNA-Seq-Workflow
Workflow for estimating gene expression in osteosarcoma RNA sequencing (RNA-Seq) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the STAR Aligner.  The version used was 2.7.6a
- Trimmomatic version 0.39
- RSEM version 1.3.3
- GATK version 4.0.12.0 (Picard included)
- All Python scripts were run on Python version 2.7.13 unless otherwise noted.  

# Workflow
**Notes:**
- Each sample has it's own directory for output files.  Each individual directory is labeled by the "Sample ID"

**Step 1) Download the hg38 genome files:** The following files were downloaded from Ensembl on 12/9/2020 and were saved to a directory named "GRCh38.p13"
- Homo_sapiens.GRCh38.dna.toplevel.fa.gz (ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/)
- Homo_sapiens.GRCh38.102.gtf.gz (ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/)

The Toplevel .fa file (Homo_sapiens.GRCh38.dna.toplevel.fa.gz) contains all sequence regions flagged as toplevel in an Ensembl schema. This includes chromsomes, regions not assembled into chromosomes and N padded haplotype/patch regions.

The .gtf file (Homo_sapiens.GRCh38.102.gtf.gz): The associated README file does not explain all the different .gtf options so I chose the one that seemed consistent with the .fa file

Both files were then unzipped as follows:
```
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.102.gtf.gz
```
**Step 2) Generate the STAR genome index files:** STAR run with the "GenomeGenerate" mode will take the user supplied reference geneome sequences (FASTA files) and annotations (GTF file) and generate genome indexes that are utilized in the mapping step.  The genome indexes are saved to disk and need only be generated once for each genome/annotation combination.
```
GENOME_DIR=<path to directory where you want the STAR index files to be saved>
FASTA=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens.GRCh38.dna.toplevel.fa"
GTF=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens.GRCh38.102.gtf"

STAR \
    --runThreadN 32 \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM 170000000000 \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles $FASTA \
    --sjdbGTFfile $GTF \
```
**Step 3) Generate the RSEM index files:** RSEM (RNA-Seq by Expectation Maximization) is a program for estimating gene and isoform expression levels from RNA-Seq data.  See https://github.com/deweylab/RSEM for details.

RSEM can extract reference transcripts from a genome if you provide it with gene annotations in a GTF/GFF file.  Alternatively, you can provide RSEM with transcript sequences directly.

The reference_fasta_file(s) is provided as either a comma-separated list of Multi-FASTA formatted files OR a directory name.  If a directory name is specified, RSEM will read all files with suffix .fa or .fasta in this directory. The reference name is the name of the reference used (in this case "human_ensembl").  RSEM will generate several reference-related files that are prefixed by this name.  This name can contain path info (e.g., "/ref/human_ensembl").

The "--gtf" option tells RSEM to assum that "reference_fasta_file(s)" contains the sequence of a genome, and it will extract transcript reference sequences using the gene annotations specified in the GTF file.
```
RSEM_INDEX_DIR=<path to directory where you want the RSEM index files to be saved>
FASTA=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens.GRCh38.dna.toplevel.fa"
GTF=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens.GRCh38.102.gtf"

rsem-prepare-reference \
    --gtf $GTF \
    $FASTA \
    $RSEM_INDEX_DIR/human_ensembl
```
**Step 4) Revert .bam files:** The downloaded osteo RNA-Seq files have already been processed and aligned. In order to use your own analysis workflow, the files first must be reverted back to their pre-processed and unmapped form. To accomplish this, we first use the Picard function, RevertSam, to take an aligned .bam file and remove alignment information in order to generate an unmapped BAM (uBAM). Details on this part of the workflow can be found in the GATK Tutorial #6484: (How to) Generate an unmapped BAM from FASTQ or aligned BAM. 
https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
It removes alignment information and any recalibrated base quality information.  It makes it possible to re-analyze the file using your own pipeline.

Standard tags cleared by default are NM, UQ, PG, MD, MQ, SA, MC and AS.  Additionally, the OQ tag is removed by the default "RESTORE_ORIGINAL_QUALITIES" parameter.  Any nonstandard tags should be removed.  To list all tags within a BAM, use the command: `samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}'` You should leave the RG tag.

The "SANITIZE" option removes reads that cause problems for certain tools such as MarkIlluminaAdapters. It removeds paired reads with missing mates, duplicate records and records with mismatches in length of bases and qualities.

For paired read files, because each read in a pair has the same query name, they are listed consecutively. We make sure to alter the previous sort order.  Coordinate sorted reads result in the aligner incorrectly estimating insert size from blocks of paired reads as they are not randomly distributed.

I use the "OUTPUT_BY_READGROUP=true" option in order to create a separate file for each readgroup.  This is neccessary because there are two different MarkDuplicate steps (one before merging each readgroup and one after merging).  This ensures that both optical and PCR duplicates are identified (see GATK tutorial 6483) The "MAX_DISCARD_FRACTION" option is informational only.  It does nto affect processing. The "SORT_ORDER=queryname", "RESTORE_ORIGINAL_QUALITIES=true", REMOVE_DUPLICATE_INFORMATION=true" and "REMOVE_ ALIGNMENT_INFORMATION=true" options are all default but I kept them in the code anyway.
```
ALIGNMENT_RUN=<Sample ID>
INPUT=<path to inptu directory>"/"$ALIGNMENT_RUN".RNA-Seq.bam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar RevertSam \
    I=$INPUT \
    O=$OUTPUT_DIR \
    OUTPUT_BY_READGROUP=true \
    SANITIZE=true \
    MAX_DISCARD_FRACTION=0.005 \
    ATTRIBUTE_TO_CLEAR=XS \
    ATTRIBUTE_TO_CLEAR=XA \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=true \
    REMOVE_DUPLICATE_INFORMATION=true \
    REMOVE_ALIGNMENT_INFORMATION=true \
    TMP_DIR=<path to temp directory>/working_temp_rs
```
**Step 5) Mark adapters:** The Picard function, MarkIllumiaAdapters, is used to take an uBAM file and rewrite it with new adapter-trimming tags.  Per tutorial 6483 on the GATK website: https://software.broadinstitute.org/gatk/documentation/article?id=6483 This is needed so that the sequenced adapters contribute minimally to the alignments.  The tool adds an "XT" tag to a read record to mark the 5' start position of the specified adapter sequence.  It also produces a metrics file.
```
#for n lanes
#The input files are the .bam files created in Step 4.  There are the same number of input .bam files as there are lanes.  

ALIGNMENT_RUN=<Sample ID>
INPUT=<path to input directory>"/"$ALIGNMENT_RUN"/"<name of .bam file created in step 4>
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar MarkIlluminaAdapters \
    I=$INPUT \
    O=$OUTPUT_DIR/uBAM_markedAdapters_lane_<n>.bam \
    M=$OUTPUT_DIR/adapter_metrics_lane_<n>.txt \
    TMP_DIR=<path to temp directory>/working_temp_miat
```
**Step 6) Convert uBAM to FASTQ:** The Picard function, SamToFastq, is used to take the uBAM files which have been modifed so that Illumina adapter sequences are marked with the XT tag and converts them to fastq files for further processing.  It produces a .fastq file in which all extant meta data (read group info, alignment info, flags and tags) are purged.  What remains are the read query names prefaced with the @ symbol, read sequences and read base quality scores.  The meta data will be added back later in the MergeBam step.  GATK actually pipes this along with the BWA alignment step and the MergeBam step but I'm still working on piping with Slurm. See Tutorial 6483 for details: 
https://software.broadinstitute.org/gatk/documentation/article?id=6483

By setting the CLIPPING_ATTRIBUTE to "XT" and by setting the CLIPPING_ACTION to 2, the program effectively removes the previously marked adapter sequences by changing their quality scores to two.  This makes it so they don't contribute to downstream read alignment or alignment scoring metrics.

The NON_PF=true option is set per the GATK tutorial.  This tells the program to retain reads marked with the 0x200 flag bit that denotes reads that do not pass quality controls (reads failing platform or vendor quality checks).
```
#For n lanes

ALIGNMENT_RUN=<Sample ID>
INPUT=<path to input directory>"/"$ALIGNMENT_RUN"/uBAM_markedAdapters_lane_<n>.bam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar SamToFastq \
    I=$INPUT \
    FASTQ=$OUTPUT_DIR/R1_lane_<n>.fastq \
    SECOND_END_FASTQ=$OUTPUT_DIR/R2_lane_<n>.fastq \
    CLIPPING_ATTRIBUTE=XT \
    CLIPPING_ACTION=2 \
    NON_PF=true \
    TMP_DIR=<path to temp directory>/working_temp_stf
```
**Step 7) Trim adapter sequences prior to alignment:** Trimmomatic is a flexibel read trimming tool for Illumina NGS data.  It performs a variety of useful trimming tasks for illumina paired-end and single ended data. See: http://www.usadellab.org/cms/?page=trimmomatic

The "PE" tells trimmomatic to run in Paired End Mode

The "trimlog file" creates a log of all read trimmings, indicating the following details:
- The read name
- The surviving sequence length
- The location of the first surviving base, aka, the amount trimmed from the start
- The location of the last surviving base in the original read
- The amount trimmed from the end

The "-phred33" option tells Trimmomatic to use phred33 instead of phred64.  Most newer sequencing data is phred33.  If you see symbols like #, !, $, % or any numbers, you know you've got phred33.

For paired-end data, two input filres are specified, and 4 output files, 2 for the 'paired' output where both reads survivied the processing and 2 for corresponding 'unpaired' output where a read survivied, but the partner read did not.

The "ILLUMINACLIP" command cuts adapter and other illumina sequences from the read. The numbers following "ILLUMINACLIP" are as followes - :2:30:12, which means 2 seed mismatches:30 palindrome clip threshold: 12 simple clip threshold.
```
#For n lanes

ALIGNMENT_RUN=<Sample ID>
ADAPTERS=<path to Trimmomatic adapters directory>"/TruSeq3-PE.fa"
WORKING_DIR=<path to input directory>"/"$ALIGNMENT_RUN
TRIMLOG=$WORKING_DIR"/Trimlog.txt"

INPUT1=$WORKING_DIR"/R1_lane_<n>.fastq"
INPUT2=$WORKING_DIR"/R2_lane_<n>.fastq"
OUTPUT1=$WORKING_DIR"/R1_0<n>.IlluminaAdapterTrimming.fastq.gz"
OUTPUT2=$WORKING_DIR"/R1_0<n>.IlluminaAdapterTrimming.unpaired1.fastq.gz"
OUTPUT3=$WORKING_DIR"/R2_0<n>.IlluminaAdapterTrimming.fastq.gz"
OUTPUT4=$WORKING_DIR"/R2_0<n>.IlluminaAdapterTrimming.unpaired1.fastq.gz"

java -Xmx48g -Xms48g -Xss2m -jar trimmomatic-0.39.jar PE \
    -trimlog $TRIMLOG \
    -phred33 \
    -threads 7 \
    $INPUT1 \
    $INPUT2 \
    $OUTPUT1 \
    $OUTPUT2 \
    $OUTPUT3 \
    $OUTPUT4 \
    ILLUMINACLIP:$ADAPTERS:2:30:12
```
**Step 8) Combine lanes prior to STAR alignment:**
```
ALIGNMENT_RUN=<Sample ID>
WORKING_DIR=<path to input directory>"/"$ALIGNMENT_RUN

cat $WORKING_DIR/R1_*.IlluminaAdapterTrimming.fastq.gz > $WORKING_DIR/R1_all_lanes.fastq.gz
cat $WORKING_DIR/R2_*.IlluminaAdapterTrimming.fastq.gz > $WORKING_DIR/R2_all_lanes.fastq.gz
```
**Step 9) Align reads with STAR:** 
The "quantMode GeneCounts" option will count the number of reads per gene while mapping. A read is counted if it overlaps (1nt or more one and only one gene.  Both ends of the paired-end read are checked for overlaps.  This requires annotations (GTF with -sjdbGTFfile option) used at the genome generation step or at the mapping step.  STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options
- column 1: gene ID
- column 2: counts for unstranded RNA-seq
- column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
- column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

The "--quantMode TranscriptomeSAM" options tells STAR to output alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam file (in addition to alignmgnets in genomic coordinates in Aligned.sam/bam files).  These transcriptomic alignments can be used with various transcript quantification software that require reads to be mapped to transcriptome, such as RSEM or eXpress.  For example, RSEM command line would look like this:
```
    rsem-calculate-expression . . . --bam Aligned.toTranscriptome.out.bam
    /path/to/RSEM/reference RSEM
```
Note, STAR first aligns reads to the entire genome, and only then searches for concordance between alignments and transcripts.  This approach offers certain advantages compared to the alignment to transcriptome only by not forcing the alignments to annotate transcripts. Using "--quantMode TranscriptomeSAM GeneCounts" will output both the Aligned.toTranscriptome.out.bam and the ReadsPreGene.out.tab outputs.
```
ALIGNMENT_RUN=<Sample ID>
INDEX=<path to directory containing the STAR index files created in Step 2>
INFILE1=<path to input directory>"/"$ALIGNMENT_RUN"/R1_all_lanes.fastq.gz"
INFILE2=<path to input directory>"/"$ALIGNMENT_RUN"/R2_all_lanes.fastq.gz"
OUT_FILE=<path to output directory>"/"$ALIGNMENT_RUN"/STAR_aligned"

STAR \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 6 \
    --outSAMunmapped None \
    --genomeDir $INDEX \
    --twopassMode Basic \
    --outFileNamePrefix $OUT_FILE \
    --outReadsUnmapped Fastq \
    --readFilesIn $INFILE1 $INFILE2 \
    --readFilesCommand zcat \
    --outSAMattributes All \
    --outFilterMultimapNmax 10 \
    --quantMode TranscriptomeSAM GeneCounts \
    --runThreadN 12
```
**Step 10) Estimate gene expression using RSEM:** 
RSEM (RNA-Seq by Expectation Maximization) is a program for estimating gene and isoform expression levels from RNA-Seq data.  See https://github.com/deweylab/RSEM for details.

Note: This step is not actually neccessary for the differential expression analysis carried out in Steps 11-12.

RSEM can extract reference transcripts from a genome if you provide it with gene annotations in a GTF/GFF3 file.  Alternatively, you can provide RSEM with transcript sequences directly.

Format: `rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name`

The "input" is the SAM/BAM/CRAM formatted input file.  RSEM requires all alignments of the same read group together.  For paired-end reads, RSEM also requires the two mates of any alignmnet be adjacent. In addition, RSEM does not allow the SEQ and QUAL fields to be empty.

The "reference_name" is the name of the reference used.  This is the output from "rsem-prepare-reference".

The "sample_name" is the name of the sample being analyzed.  All output files are prefixed by this name.
```
SAMPLE=<sample name>
RSEM_INDEX_DIR=<path to directory containing the RSEM index files created in Step 3>
INPUT_FILE=<path to input directory>"/"$SAMPLE"/STAR_alignedAligned.toTranscriptome.out.bam"
OUTPUT_PREFIX=<path to output directory>"/"$SAMPLE"/RSEM_expression"

rsem-calculate-expression \
    --paired-end \
    --alignments \
    -p 32 \
    $INPUT_FILE \
    $RSEM_INDEX_DIR/human_ensembl \
    $OUTPUT_PREFIX
```
**Step 11) Prepare gene count data to be inputed into DESeq2:**  The data from the STAR alignment step will be prepared for input into DESeq2 using the Python program, DESeq2_input.py.  
**Step 12) Identify differentially expressed genes between samples with and without nuclear envelope or nuclear pore complex pathway aberrations:**  The R program, DESeq2_Analysis.R, is used to identify differentially expressed genes using DESeq2.  
