import gzip
import shutil
import argparse
import tempfile
import pyBigWig
from pyfaidx import Fasta

def reverse_complement(sequence):
    complement_dict = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'n':'n', 'N':'N'}
    complemented_sequence = [complement_dict[base] for base in reversed(sequence)]
    reverse_complement = ''.join(complemented_sequence)
    return reverse_complement

def filter_mappability(final_output_file, bw, min_mappability_score):
    lines_to_drop = []
    with open(final_output_file, "r") as final_output:
        for i, line in enumerate(final_output):
            fields = line.strip().split("\t")
            chrom = fields[0]
            end = int(fields[1])
            mappability = bw.stats(chrom, end-5, end+5)[0]
            if mappability is not None:
                if mappability < min_mappability_score:
                    lines_to_drop.append(i)
        print(lines_to_drop)
    lines_to_drop_set = set(lines_to_drop)
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file,  open(final_output_file, "r") as final_output:
        for i, line in enumerate(final_output):
                if i not in lines_to_drop_set:
                    temp_file.write(line)
    shutil.move(temp_file.name, final_output_file)

def blacklist_bed(final_output_file, blacklist_file):
    lines_to_drop = []
    with gzip.open(blacklist_file, "rt") as f,  open(final_output_file, "r") as final_output:
        for blacklist_line in f:
            values = blacklist_line.strip().split("\t")
            chrom_blacklist = values[0]
            start_range = int(values[1])
            end_range = int(values[2])
            for i, line in enumerate(final_output):
                fields = line.strip().split("\t")
                chrom = fields[0]
                end = int(fields[1])
                if chrom==chrom_blacklist:
                        if (start_range <= end-5 <= end_range 
                            or start_range <= end+5 <= end_range
                            or end-5 <= start_range <= end+5
                            or end-5 <= end_range <= end+5
                        ):
                            lines_to_drop.append(i)

    lines_to_drop_set = set(lines_to_drop)
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file,  open(final_output_file, "r") as final_output:
        for i, line in enumerate(final_output):
                if i not in lines_to_drop_set:
                    temp_file.write(line)
    shutil.move(temp_file.name, final_output_file)

def get_sequence(final_output_file, fasta):
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file, open(final_output_file, "r") as final_output:
        for line in final_output:
            fields = line.strip().split("\t")
            chromosome = fields[0]
            end = int(fields[2])
            methylation = str(fields[6])+","+str(fields[7])
            sequence_pos = str(fasta[str(chromosome)][(end-6):(end+5)])
            sequence_neg = str(fasta[str(chromosome)][(end-5):(end+6)])
            sequence_neg = reverse_complement(sequence_neg)
            temp_file.write(f"{chromosome}\t{end}\t{methylation}\t{sequence_pos}\t{sequence_neg}\n")
    shutil.move(temp_file.name, final_output_file)

def remove_duplicates(source_bed_file, final_output_file):
    with open(source_bed_file, "r") as source_bed:
        lines = source_bed.readlines()
        num_lines = len(lines)
        i = 1
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            while i < num_lines:
                line = lines[i].strip().split("\t")
                if i + 1 < num_lines:
                    next_line = lines[i + 1].strip().split("\t")
                    if line[2] == next_line[1]:
                        methylation = (float(next_line[6])/100 * int(next_line[7]) + float(line[6])/100 * int(line[7])) / (float(int(line[7]) + int(next_line[7])))
                        line[6] = str(methylation*100)
                        line[7] = str(int(line[7])+int(next_line[7]))
                        i += 1
                if line[5]=="-":
                    line[2]=str(int(line[2])-1)
                    line[5]="+"
                temp_file.write("\t".join(line) + "\n")
                i += 1
    shutil.move(temp_file.name, final_output_file)

def isolate_c_g(source_bed_file, blacklist_file, mappability_file, min_mappability_score, reference_genome_file, final_output_file):
    fasta = Fasta(reference_genome_file)
    bw = pyBigWig.open(mappability_file)
    remove_duplicates(source_bed_file, final_output_file)
    get_sequence(final_output_file, fasta)
    blacklist_bed(final_output_file, blacklist_file)
    filter_mappability(final_output_file, bw, min_mappability_score)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source_bed_file', type=str) ## The source .bed file containing the locations of CG sites.
    p.add_argument('--blacklist_file', type=str) ## The .bed file containing locations to blacklist.
    p.add_argument('--mappability_file', type=str) ## A .bw file containing a mappability track
    p.add_argument("--min_mappability_score", type=float) # Minimum mappability score to include.
    p.add_argument('--reference_genome_file', type=str) ## A .fa file containing the reference human genome.
    p.add_argument('--final_output_file', type=str) ## The path and name of the final output file.
    args = p.parse_args()
    isolate_c_g(**vars(args))




