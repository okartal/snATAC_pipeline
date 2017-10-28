Single-cell / single nuclei ATAC-seq  pipeline 
===================================================
!!! For the old version of the shell scripts, check the legend README [here](./README_legend.md). 

Below is the use of this pipeline wrappered using bigdatascript (bds) language.

# Installation
## For tscc user

The dependent softwares and libs are installed in a conda enrionment that can directly be loaded. 
Add the following into your `~/.bashrc`
```bash
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH
export PATH="$PATH:/projects/ps-epigen/software/miniconda3/bin/"
export _JAVA_OPTIONS="-Xms256M -Xmx728M -XX:ParallelGCThreads=1"
export PATH="$PATH:/projects/ps-epigen/software/.bds/"
export PICARDROOT="/projects/ps-epigen/software/miniconda3/envs/bds_atac/share/picard-1.126-4/"
```
## For other user 
1. Install conda: check https://github.com/kundajelab/atac_dnase_pipelines#conda
2. install bds: check https://github.com/kundajelab/atac_dnase_pipelines#bigdatascript
3. Run `bash ./install_dependencies.sh` to generate `bds_scATAC` environment to encapsulate dependent softwares. 
    
# Usage 
``` shell
$ scATAC.bds -help
Picked up _JAVA_OPTIONS: -Xms256M -Xmx728M -XX:ParallelGCThreads=1
== scATAC pipeline settings
	-barcode_dir <string>        : folder that contains r7_ATAC, i7_ATAC, i5_ATAC and r5_ATAC barcode.
	-bowtie2_idx <string>        : Bowtie2 indexed reference genome.
	-mark_duplicate <string>     : path to picard MarkDuplicates.jar [MarkDuplicates.jar].
	-max_barcode_mismatch <int>  : max barcode mismatch allowed for barcode error correction [2].
	-min_read <int>              : cells with reads less than 'min_read' will be filtered [500].
	-prefix <string>             : prefix of output files.
	-r1 <string>                 : fastq.gz file that contains forward reads (only .gz allowed).
	-r2 <string>                 : fastq.gz file that contains reverse reads (only .gz allowed).
	-threads <int>               : Set threads for Bowtie2 mapping, sorting and duplication removal [1].
```

## For tscc user 
Check the `pbs` files in the `\examples` folder. 

## Input 
Currently it run start from the fastq files after decomplex and demultiplex. 

# Licence
MIT






