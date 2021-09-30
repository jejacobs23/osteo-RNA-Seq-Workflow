# osteo-RNA-Seq-Workflow
Workflow for estimating gene expression in osteosarcoma RNA sequencing (RNA-Seq) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the STAR Aligner.  The version used was 
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
