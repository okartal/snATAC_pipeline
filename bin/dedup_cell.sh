#!/bin/bash
PREFIX=$1; THREADS=$2; 

iterator=0
for line in `ls $PREFIX\_tmp/cells | grep .sam`
do
	barcode="${line%.*}"
	samtools view -bS $PREFIX\_tmp/cells/$barcode.sam \
	| samtools sort -@ $THREADS - -o $PREFIX\_tmp/cells/$barcode.sorted.bam
	samtools rmdup $PREFIX\_tmp/cells/$barcode.sorted.bam $PREFIX\_tmp/cells/$barcode.sorted.nodup.bam
    iterator=$((iterator + 1))
    if (( $iterator % $THREADS == 0 ))           # no need for brackets
    then
  	  while kill -0 $! 2>/dev/null
  	  do
  	    sleep 5
  	  done
    fi
done
