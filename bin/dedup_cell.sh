#!/bin/bash
PREFIX=$1; THREADS=$2; 

iterator=0
ls $PREFIX\_tmp/cells | grep .sam | sed "s/\.sam//g" | xargs -i --max-procs=$THREADS \
                                          bash -c "samtools view -bS $PREFIX\_tmp/cells/{}.sam \
	    | samtools sort - -o $PREFIX\_tmp/cells/{}.sorted.bam  >/dev/null 2>&1; \
samtools rmdup $PREFIX\_tmp/cells/{}.sorted.bam $PREFIX\_tmp/cells/{}.sorted.nodup.bam >/dev/null 2>&1"


