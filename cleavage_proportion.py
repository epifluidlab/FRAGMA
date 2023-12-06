import pysam
import shutil
import argparse
import tempfile
import concurrent.futures

CHUNK_SIZE = 3000 

def process_chunk(chunk, tbx_chunk, minimum_fragments):
    processed_lines = []
    for line in chunk:
        result = process_line(line, tbx_chunk, minimum_fragments)
        if result is not None:
            processed_lines.append(result)
    return processed_lines

def fetch_fragment_information(tbx_fragment, chrom, location, positive):
    fragments=[]
    tbx_fragment = pysam.TabixFile(tbx_fragment)
    for row in tbx_fragment.fetch(chrom, location-6, location+6):
        fields = str(row).split('\t')
        if positive:
            if fields[5]=="+":
                fragments.append(str(row))
        else:
            if fields[5]=="-":
                fragments.append(str(row))
    num_fragments=len(fragments)
    return fragments, num_fragments

def process_line(line, tbx_line, minimum_fragments):
    fields = line.strip().split("\t")
    chrom = int(fields[0])
    location_pos = int(fields[1])
    location_neg = location_pos + 1
    methylation_line = fields[2].split(',')
    total_line = float(methylation_line[1])
    methylation = float(methylation_line[0])
    pos_strand = fields[3]
    neg_strand = fields[4]
    negative_array = [None] * 11
    positive_array = [None] * 11
    for i in range(-5,6):
        positive_array[i]=0.0
        negative_array[i]=0.0
    if methylation >= 70 or methylation <= 30:
        pos_fragments, total_fragments_positive = fetch_fragment_information(tbx_line, chrom, location_pos, True)
        neg_fragments, total_fragments_negative = fetch_fragment_information(tbx_line, chrom, location_neg, False)
        if total_fragments_positive > minimum_fragments and total_fragments_negative > minimum_fragments:
            for row in pos_fragments:
                start_end = row.strip().split("\t")
                start = int(start_end[1])
                end = int(start_end[2])
                if location_pos - 5 <= start <= location_pos + 5:
                    relative_location = start - location_pos
                    positive_array[relative_location] += 1
                if location_pos - 5 <= end <= location_pos + 5:
                    relative_location = end - location_pos
                    positive_array[relative_location] += 1

            for row in neg_fragments:
                start_end = row.strip().split("\t")
                start = int(start_end[1])
                end = int(start_end[2])
                if location_neg - 5 <= start <= location_neg + 5:
                    relative_location = start - location_neg
                    negative_array[relative_location] += 1
                if location_neg - 5 <= end <= location_neg + 5:
                    relative_location = end - location_neg
                    negative_array[relative_location] += 1

            positive_array_ordered=[]
            negative_array_ordered=[]
            for i in range(-5, 6):
                positive_array_ordered.append(positive_array[i])
                negative_array_ordered.append(negative_array[i])
            return f"{chrom}\t{location_pos}\t{methylation}\t{total_line}\t{pos_strand}\t{neg_strand}\t{positive_array_ordered}\t{total_fragments_positive}\t{negative_array_ordered}\t{total_fragments_negative}\n"

    return None
    
def cleavage_proportion(cg_file, frag_file, minimum_fragments):
    tbx = frag_file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
        with open(cg_file, "r") as cg_sites:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = []
                chunk = []

                for i, line in enumerate(cg_sites):
                    chunk.append(line)
                    if i % CHUNK_SIZE == 0:
                        futures.append(executor.submit(process_chunk, chunk.copy(), tbx, minimum_fragments))
                        chunk = []

                if chunk:
                    futures.append(executor.submit(process_chunk, chunk.copy(), tbx, minimum_fragments))

                for future in concurrent.futures.as_completed(futures):
                    results = future.result()
                    for result in results:
                        temp_file.write(result)

        shutil.move(temp_file.name, cg_file)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--cg_file', type=str) ## The source .bed file containing the locations of CG sites.
    p.add_argument('--frag_file', type=str) ## The .bed file containing fragments.
    p.add_argument('--minimum_fragments', type=int, default=0) ## The minimum fragments to include CG site in cleavage proportion list.
    args = p.parse_args()
    cleavage_proportion(**vars(args))



