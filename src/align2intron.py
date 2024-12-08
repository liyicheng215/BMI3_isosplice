import pysam
import numpy as np


# check the matching pattern near the intron
def find_matching_pattern(record, n_pos, k):
    cigar = record.cigartuples
    # check the matching pattern
    mismatches = ""
    counts = 0
    for i in n_pos:
        count = cigar[i][1]
        if count >= k-counts:
            mismatches += "".join([str(cigar[i][0])] * (k-counts))
            break
        else:
            counts += count
            mismatches += "".join([str(cigar[i][0])] * cigar[i][1])
    return mismatches


# score the introns based on the mismatch near the intron, the higher the score, the more mismatch/nearer, less likely to be SJ
def score_introns(record, n_pos, k):
    mismatches = find_matching_pattern(record, n_pos, k)
    mismatch_score = 0
    for i in range(len(mismatches)):
        if mismatches[i] != "0":
            mismatch_score += k-i
    return mismatch_score


# find all introns in the alignment record, return their start and end positions on the reference genome and score
def find_introns(record, k):
    introns = []
    crm = record.reference_name  # chromatin
    if record.is_reverse:
        strand = "-"
    else:
        strand = "+"
    gene_pos = record.reference_start  # use the start position of the reads on the chromosome as TSS
    ref_pos = record.reference_start  # the position of the reads on the current reference sequence, used to find the start and end positions of the intron

    numN = record.cigarstring.count("N")  # how many Ns in a read?
    if numN == 1:
        onlyN = 1
    else:
        onlyN = 0

    cigar = record.cigartuples
    for i in range(len(cigar)):
        if cigar[i][0] == 3:  # 3=N=the bases spanned, representing introns
            start = ref_pos
            end = ref_pos + cigar[i][1] - 1
            mismatch_score = (score_introns(record, range(i-1, -1, -1), k) +
                              score_introns(record, range(i+1, len(cigar), 1), k))   # mismatch situation and score before and after the intron
            introns.append((crm, strand, gene_pos, start, end, onlyN, mismatch_score))
        if cigar[i][0] in [0, 2, 3]:  # 0=M=matching bases, 2=D=missing bases from the reference genome, both will move the position of the reference sequence
            ref_pos += cigar[i][1]
    return introns


# read the bam file and keep only the best alignment
def filter_best_alignment(input_bam, k):
    alignments = {}
    max_qualities = {}
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        for record in infile:
            flag = record.flag
            if not (flag & 0x904):   # exclude unaligned, secondary, PCR or optical duplicates and supplementary alignments
                qname = record.query_name
                mapq = record.mapping_quality
                if qname not in alignments.keys() or mapq > max_qualities[qname]:
                    max_qualities[qname] = mapq
                    intron = find_introns(record, k)
                    if intron:
                        alignments[qname] = intron
    return alignments


# convert the alignment to introns
def align2intron(input_bam, k):
    introns = []
    alignment = filter_best_alignment(input_bam, k)
    index = 0
    read_id = 0
    for read, intron in alignment.items():
        read_id += 1
        intron_id = 0
        for itr in intron:
            index += 1
            intron_id += 1
            new_intron = ["index"+str(index)] + ["read"+str(read_id)] + ["intron"+str(intron_id)] + list(itr)
            introns.append(new_intron)
    return introns


if __name__ == '__main__':
    Introns = align2intron('../../align2intron/ce11_ONT_20_sorted.bam', k=8)
    np.savetxt('1206Introns.csv', Introns, delimiter=",", fmt='%s')
    print(Introns)
