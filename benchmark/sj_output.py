from sys import stderr, exit
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType
import csv

def print_stats(exon_lengths, intron_lengths, trans_lengths, genes, trans, file=stderr):
    print('genes: {}, genes with multiple isoforms: {}'.format(
        len(genes), sum(len(v) > 1 for v in genes.values())),
        file=file)
    print('transcripts: {}, transcript avg. length: {:.0f}'.format(
        len(trans), sum(trans_lengths.elements()) // len(trans)),
        file=file)
    print('exons: {}, exon avg. length: {:.0f}'.format(
        sum(exon_lengths.values()),
        sum(exon_lengths.elements()) // sum(exon_lengths.values())),
        file=file)
    print('introns: {}, intron avg. length: {:.0f}'.format(
        sum(intron_lengths.values()),
        sum(intron_lengths.elements()) // sum(intron_lengths.values())),
        file=file)
    print('average number of exons per transcript: {:.0f}'.format(
        sum(exon_lengths.values()) // len(trans)),
        file=file)

def extract_splice_sites(gtf_file, verbose=False, csv_file=None):
    genes = dd(list)
    trans = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and calculate introns
    introns = []
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        valid_introns = []
        for i in range(1, len(exons)):
            intron_start = exons[i - 1][1] + 1
            intron_end = exons[i][0] - 1
            intron_length = intron_end - intron_start + 1
            if intron_length >= 4:  # Remove introns shorter than 4 bp
                valid_introns.append((intron_start, intron_end, strand))
        introns.extend([(chrom, start, end, strand) for start, end, strand in valid_introns])

    # Print the unique splice sites (junctions)
    for chrom, start, end, strand in sorted(introns):
        # Zero-based offset
        print(f'{chrom}\t{start - 1}\t{end - 1}\t{strand}')

    # Print some stats if asked or if csv_file is provided
    if verbose or csv_file:
        exon_lengths, intron_lengths, trans_lengths = \
            Counter(), Counter(), Counter()
        for chrom, strand, exons in trans.values():
            tran_len = 0
            for i, exon in enumerate(exons):
                exon_len = exon[1] - exon[0] + 1
                exon_lengths[exon_len] += 1
                tran_len += exon_len
                if i == 0:
                    continue
                intron_len = exons[i][0] - exons[i - 1][1] - 1
                if intron_len >= 4:  # Include only introns >= 4 bp
                    intron_lengths[intron_len] += 1
            trans_lengths[tran_len] += 1

        if csv_file:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Statistic', 'Value'])
                writer.writerow(['Genes', len(genes)])
                writer.writerow(['Genes with multiple isoforms', sum(len(v) > 1 for v in genes.values())])
                writer.writerow(['Transcripts', len(trans)])
                writer.writerow(['Transcript avg. length', sum(trans_lengths.elements()) // len(trans)])
                writer.writerow(['Exons', sum(exon_lengths.values())])
                writer.writerow(['Exon avg. length', sum(exon_lengths.elements()) // sum(exon_lengths.values())])
                writer.writerow(['Introns', sum(intron_lengths.values())])
                writer.writerow(['Intron avg. length', sum(intron_lengths.elements()) // sum(intron_lengths.values())])
                writer.writerow(['Average number of exons per transcript', sum(exon_lengths.values()) // len(trans)])

                # 写入内含子详细信息
                writer.writerow([])
                writer.writerow(['Chromosome', 'Start', 'End', 'Strand'])

                seen_introns = set()
                for chrom, start, end, strand in sorted(introns):
                    intron = (chrom, start, end, strand)
                    if intron not in seen_introns:
                        writer.writerow(intron)
                        seen_introns.add(intron)


        else:
            print_stats(exon_lengths, intron_lengths, trans_lengths, genes, trans)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract splice junctions from a GTF file')
    parser.add_argument('gtf_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GTF file (use "-" for stdin)')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('--csv',
                        default='splice_sites_statistics.csv',  # Setting default values
                        help='output CSV file for statistics')

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    extract_splice_sites(args.gtf_file, args.verbose, args.csv)