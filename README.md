# osteo-RNA-Seq-Workflow
Workflow for estimating gene expression in osteosarcoma RNA sequencing (RNA-Seq) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the STAR Aligner.  The version used was 2.7.6a
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
**Step 3) :** The Picard function, RevertSam, is ued to take an aligned .bam file and remove alignment information in order to generate an unmapped BAM (uBAM). Details on this part of the workflow can be found in the GATK Tutorial #6484: (How to) Generate an unmapped BAM from FASTQ or aligned BAM. 
https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
It removes alignment information and any recalibrated base quality information.  It makes it possible to re-analyze the file using your own pipeline.

Standard tags cleared by default are NM, UQ, PG, MD, MQ, SA, MC and AS.  Additionally, the OQ tag is removed by the default "RESTORE_ORIGINAL_QUALITIES" parameter.  Any nonstandard tags should be removed.  To list all tags within a BAM, use the command: samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}' You should leave the RG tag.

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
**Step 4) :** The Picard function, MarkIllumiaAdapters, is used to take an uBAM file and rewite it with new adapter-trimming tags.  Per tutorial 6483 on the GATK website: https://software.broadinstitute.org/gatk/documentation/article?id=6483 This is needed so that the sequenced adapters contribute minimally to the alignments.  The tool adds an "XT" tag to a read record to mark the 5' start position of the specified adapter sequence.  It also produces a metrics file.
```
#for n lanes
#The input files are the .bam files created in Step .  There are the same number of input .bam files as there are lanes.  

ALIGNMENT_RUN=<Sample ID>
INPUT=<path to input directory>"/"$ALIGNMENT_RUN"/"<name of .bam file created in step >
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar $PICARD_DIR/picard.jar MarkIlluminaAdapters \
    I=$INPUT \
    O=$OUTPUT_DIR/uBAM_markedAdapters_lane_<n>.bam \
    M=$OUTPUT_DIR/adapter_metrics_lane_<n>.txt \
    TMP_DIR=<path to temp directory>/working_temp_miat
```


**Step 5) :** The Picard function, SamToFastq, is used to take the uBAM files which have been modifed so that Illumina adapter sequences are marked with the XT tag and converts them to fastq files for further processing.  It produces a .fastq file in which all extant meta data (read group info, alignment info, flags and tags) are purged.  What remains are the read query names prefaced with the @ symbol, read sequences and read base quality scores.  The meta data will be added back later in the MergeBam step.  GATK actually pipes this along with the BWA alignment step and the MergeBam step but I'm still working on piping with Slurm. See Tutorial 6483 for details: 
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
