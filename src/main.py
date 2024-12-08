import argparse
from align2intron import align2intron
from graph_simplification import simplifyDAG
from build_DAG import build_DAG
from score_SJ import train_and_score_sj

def main():
    # set up command line parsing
    parser = argparse.ArgumentParser(description="Intron and SJ finder from BAM file")
    subparsers = parser.add_subparsers(dest="command")

    # handle findSJ command
    findSJ_parser = subparsers.add_parser("findSJ", help="Find splice junctions from BAM file")
    findSJ_parser.add_argument("input_bam", help="Input BAM file")
    findSJ_parser.add_argument("--k", type=int, default=8, help="Score threshold for mismatches")

    # Parser for simplifyDAG command
    simplify_parser = subparsers.add_parser("simplifyDAG", help="Simplify DAG based on splice junctions")
    simplify_parser.add_argument("input_file", default = "DAG_edge.csv", help="Input CSV file with the splice junction data")
    simplify_parser.add_argument("--thres", type=float, default=0.01, help="Threshold for simplification")
    simplify_parser.add_argument("--thres_dominant", type=float, default=0.99, help="Threshold for dominant edges")

    # Parser for buildDAG command
    buildDAG_parser = subparsers.add_parser("buildDAG", help="Build a DAG from intron data")
    buildDAG_parser.add_argument("input_file", default = "SJ_scores.csv", help="Input CSV file with intron data")
    buildDAG_parser.add_argument("--genome_length", type=int, default=100286401, help="Length of the genome (default: 100286401)")
    buildDAG_parser.add_argument("--cluster_shift", type=int, default=8, help="Cluster shift value (default: 8)")
    buildDAG_parser.add_argument("--gene_shift", type=int, default=8, help="Gene shift value (default: 8)")

    # Parser for scoreSJ command
    scoreSJ_parser = subparsers.add_parser('scoreSJ', help='Train and score SJ data using RandomForest')
    scoreSJ_parser.add_argument('csv_file', type=str, default = "cluster_scores.csv", help='Path to the CSV file containing SJ cluster.')
    scoreSJ_parser.add_argument('csv_label', type=str, default = "cluster_label.csv", help='Path to the CSV file containing SJ label.')
    scoreSJ_parser.add_argument('--split', type=float, default=None, help='Threshold for labeling SJ scores (0-1).')

    args = parser.parse_args()

    if args.command == "findSJ":
        # call align2intron function from findSJ.py
        print(f"Processing BAM file: {args.input_bam} with k={args.k}")
        Introns = align2intron(args.input_bam, k=args.k)
        output_file = 'SJ_scores.csv'
        # save results to a CSV file
        with open(output_file, 'w') as f:
            for intron in Introns:
                f.write(",".join(map(str, intron)) + "\n")
        print(f"Splice junctions and introns have been saved to {output_file}")
    
    elif args.command == "simplifyDAG":
        print(f"Processing input file: {args.input_file} with thresholds {args.thres} and {args.thres_dominant}")
        simplifyDAG(args.input_file, thres=args.thres, thres_dominant=args.thres_dominant)

    elif args.command == "buildDAG":
        print(f"Processing input file: {args.input_file}")
        print(f"Genome length: {args.genome_length}, Cluster shift: {args.cluster_shift}, Gene shift: {args.gene_shift}")
        build_DAG(args.input_file, genome_length=args.genome_length, cluster_shift=args.cluster_shift, gene_shift=args.gene_shift)

    elif args.command == "scoreSJ":
        train_and_score_sj(args.csv_file, args.csv_label, args.split)

if __name__ == "__main__":
    main()
