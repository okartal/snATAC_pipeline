single-cell / single nuclei ATAC-seq  pipeline 
===================================================



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

## conda install (under construct) 
Run `bash ./install_dependencies.sh` to generate `bds_scATAC` environment to encapsulate dependent softwares. 
    
# Usage 

## For tscc user 
Check the `pbs` files in the `\examples` folder. 

## Input 
Currently it run start from the fastq files after decomplex and demultiplex. 

# Licence
MIT






