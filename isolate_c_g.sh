#!/bin/sh

pip install pyBigWig
python ./isolate_c_g.py \
    --source_bed_file /ocean/projects/mcb190124p/dnaase/projects/ccinference/broad_data/wgbs_bam/HD_45.WGBS.markDuplicates.cpg.filtered.sort.CG.strand.6plus2.bed \
    --blacklist_file /jet/home/rbandaru/ocean_rbandaru/FRAGMA/data/BLACKLIST_MAPPABILITY/ENCFF001TDO.bed.gz \
    --mappability_file /jet/home/rbandaru/ocean_rbandaru/FRAGMA/data/BLACKLIST_MAPPABILITY/human_g1k_v37.45mer.mappability.bw \
    --min_mappability_score 0.9 \
    --reference_genome_file /ocean/projects/mcb190124p/shared/data/genomes/hg19/human_g1k_v37.fa \
    --final_output_file /jet/home/rbandaru/ocean_rbandaru/FRAGMA/data/HD_45_CG_FINAL_UNMASKED.bed \

    