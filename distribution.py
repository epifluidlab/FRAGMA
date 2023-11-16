import pysam
import argparse
import seaborn as sns
import matplotlib.pyplot as plt


def fetch_fragment_information(tbx_loc, chrom, location):
    fragments=[]
    for row in tbx_loc.fetch(chrom, location-5, location+5):
        fragments.append(str(row))
    num_fragments=len(fragments)
    return num_fragments

def distribution(source_bed_file, mappability_file):
    tbx = pysam.TabixFile(mappability_file)
    fragment_scores = []
    lines_processed = 0
    with open(source_bed_file, "r") as final_output:
        for line in final_output:
            fields = line.strip().split("\t")
            chrom = int(fields[0])
            end = int(fields[1])
            fragments = fetch_fragment_information(tbx, chrom, end)
            if fragments is not None:
                fragment_scores.append(fragments)
            lines_processed += 1
            if lines_processed == 10000:
                break
    sns.histplot(fragment_scores, bins=100, kde=False, color='blue')
    plt.title('Distribution of Mappability Scores')
    plt.xlabel('Number of Fragments')
    plt.ylabel('Frequency')
    plt.savefig('Fragment_Number_Scores.png')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source_bed_file', type=str) ## The source .bed file containing the locations of CG sites.
    p.add_argument('--mappability_file', type=str) ## A .bw file containing a mappability track.
    args = p.parse_args()
    distribution(**vars(args))