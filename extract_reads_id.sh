########################################################################
#####     Extract sequencing ID from bam files by samtools         #####
#####                by Yufan (Harry) Zhou at 2025/10/14          ######
########################################################################

#!/bin/bash
tissue=$1
path=$2

echo $tissue
echo $path

if [[ ! -f ${path}/${tissue}.rg.dup.sort.bam ]]; then
    echo "Error: BAM file ${path}/${tissue}.rg.dup.sort.bam not found!"
    return 1
fi
if [[ ! -f ${tissue}.variants.txt ]]; then
    echo "Error: Variants file ${tissue}.variants.txt not found!"
    return 1
fi

extract_reads() {
    local tissue=$1
    local path=$2
    while read CHROM POS REF ALT; do
        samtools view ${path}/${tissue}.rg.dup.sort.bam ${CHROM}:${POS}-${POS} \
        | awk -v c=$CHROM -v p=$POS -v r=$REF -v a=$ALT '{print c "\t" p "\t" r "\t" a "\t" $0}' \
        | tr -d '\r'
    done < ${tissue}.variants.txt > ${tissue}.read_ids.txt
}

extract_reads ${tissue} ${path}
