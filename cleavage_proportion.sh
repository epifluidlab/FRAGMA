#!/bin/sh
pip install pysam
python ./cleavage_proportion.py \
    --cg_file /jet/home/rbandaru/ocean_rbandaru/FRAGMA/fractionalized_data/part_.bed \
    --frag_file /jet/home/rbandaru/ocean_rbandaru/FRAGMA/data/HD_45.WGBS.frag.bed.gz \
    --minimum_fragments 20 \

