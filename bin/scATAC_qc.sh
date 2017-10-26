#! /bin/bash

usage(){
cat << EOF
EOF
} 

while getopts ":r:l:t:c:n:d:s:g:m:" opt;
do
	case "$opt" in
	    r) R1=$OPTARG;;
	    l) LOG=$OPTARG;;            
            t) T_BAM=$OPTARG;;
            c) C_BAM=$OPTARG;;
            n) N_BAM=$OPTARG;;
            d) D_BAM=$OPTARG;;
            s) D_BAM_S=$OPTARG;;
            g) G_BAM=$OPTARG;;                        
	    m) MIN_READ=$OPTARG;;            
		\?) usage
			echo "input error"
			exit 1
			;;
	esac
done


TOTAL_READS=$(zcat $R1 | wc -l | awk '{print $1}')
UNIQ_READS=$(samtools view $T_BAM | wc -l | awk '{print $1}')
BARCODE_CORRECTED_READS=$(samtools view $C_BAM | wc -l | awk '{print $1}')
BEFORE_PICARD_READS=$(samtools view $N_BAM | wc -l | awk '{print $1}')
AFTER_PICARD_READS=$(samtools view $D_BAM | wc -l | awk '{print $1}')
TOTAL_CELLS=$(awk -v cutoff="$MIN_READ" '{if ($2 > cutoff) print }' $D_BAM_S | wc -l | awk '{print $1}')
FINAL_READS=$(samtools view $G_BAM | wc -l | awk '{print $1}')
TOTAL_BARCODE=$(samtools view $C_BAM | awk '{split($1,a,":"); print a[1]}' | sort | uniq -c | wc -l | awk '{print $1}')


echo "================================ Summary ==================================" 2>&1 | tee -a $LOG
echo "Total number of raw reads: $((TOTAL_READS/4))" 2>&1 | tee -a $LOG
echo "Uniquely mapped reads (MAPQ>=30): $UNIQ_READS" 2>&1 | tee -a $LOG
echo "Reads left after barcode correction: $BARCODE_CORRECTED_READS" 2>&1 | tee -a $LOG
echo "Uniquely barcode combinations: $TOTAL_BARCODE" 2>&1 | tee -a $LOG
echo "Estimated PCR duplication rate: $(echo "scale=2; ($BEFORE_PICARD_READS - $AFTER_PICARD_READS) / $BEFORE_PICARD_READS * 100"| bc -l)%" 2>&1 | tee -a $LOG
echo "Total number of reads left: $FINAL_READS" 2>&1 | tee -a $LOG
echo "Number of cells with more than $MIN_READ reads: $TOTAL_CELLS" 2>&1 | tee -a $LOG

sort -k2,2n $D_BAM_S \
| awk -v cutoff="$MIN_READ" '{if ($2 > cutoff) print $2}' \
| awk '
function max(x){i=0;for(val in x){if(i<=x[val]){i=x[val];}}return i;}
function min(x){i=max(x);for(val in x){if(i>x[val]){i=x[val];}}return i;}
{
     a[x++]=$1
     b[$1]++
}
END{
	print "Min number of reads for selected cells: "min(a) 
    print "Median number of reads for selected cells: "a[int((x-1)/2)] 
    print "Max number of reads for selected cells: "max(a)
}' 2>&1 | tee -a $LOG


echo "FINISHED (`date`)"| tee -a $LOG
