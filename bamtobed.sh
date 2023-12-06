#!/bin/sh
source /opt/packages/anaconda3-2022.10/etc/profile.d/conda.sh
conda init bash
module load samtools
module load bedtools

#do for HD_45 and #HD_46
samtools sort -n -o ./HD_45_WGBS_SORTED_UNMASKED.bam -@ 4 /ocean/projects/mcb190124p/dnaase/projects/ccinference/broad_data/wgbs_bam/HD_45.WGBS.markDuplicates.bam
samtools view -bh -f 3 -F 3852 -G 48 -q 30 \
  ./HD_45_WGBS_SORTED_UNMASKED.bam |
  bedtools bamtobed -bedpe -mate1 -i stdin |
  awk -F'\t' -v OFS="\t" '{if ($1!=$4) next; if ($9=="+") {s=$2;e=$6} else {s=$5;e=$3} if (e>s) print $1,s,e,$7,$8,$9}' |
  sort -k1,1V -k2,2n |
  bgzip > ./HD_45_WGBS_SORTED_UNMASKED.frag.bed.gz
tabix -p bed ./HD_45_WGBS_SORTED_UNMASKED.frag.bed.gz