import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare reference and result CSV files for splicing junctions.')
    parser.add_argument('reference', type=str, help='Reference CSV file path')
    parser.add_argument('result', type=str, help='Result CSV file path')
    return parser.parse_args()

def load_data(reference_path, result_path):
    try:
        # Load the reference file
        with open(reference_path, 'r') as file:
            lines = file.readlines()

        # Find the start line of the data section in the reference file
        start_line = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('Chromosome'):
                start_line = i
                break

        reference = pd.read_csv(reference_path, skiprows=start_line, header=0, names=['chr', 'start', 'end', 'strand'])

        # Load the result file
        result = pd.read_csv(
            result_path,  # Path
            header=0,  # First row as column names
            usecols=[1, 2, 3, 4, 12],  # Select the required columns
            names=['chr', 'strand', 'start', 'end', 'score'],  # Rename columns
        )
    except FileNotFoundError:
        print("Error: One of the input files was not found.")
        raise
    except pd.errors.EmptyDataError:
        print("Error: One of the input files is empty or improperly formatted.")
        raise

    return reference, result



def build_reference_dict(reference):
    ref_dict = {}
    for _, row in reference.iterrows():
        chr_only = row['chr']
        if chr_only not in ref_dict:
            ref_dict[chr_only] = []
        ref_dict[chr_only].append((row['start'], row['end']))
    return ref_dict


def build_result_dict(result):
    res_dict = {}
    for _, row in result.iterrows():
        chr_only = row['chr']
        start, end, score = row['start'], row['end'], row['score']

        # If the chromosome has not been added to the dictionary yet, initialize it as a set.
        if chr_only not in res_dict:
            res_dict[chr_only] = set()

        # Use a set to store tuples of `(start, end, score)` to avoid duplicate records.
        res_dict[chr_only].add((start, end, score))

    # Convert the set to a list if subsequent operations need to be consistent with the original code.
    for chr_only in res_dict:
        res_dict[chr_only] = list(res_dict[chr_only])

    return res_dict


import matplotlib.pyplot as plt
import numpy as np

def compare(reference, result):
    reference_dict = build_reference_dict(reference)
    result_dict = build_result_dict(result)

    reference_sj = sum(len(v) for v in reference_dict.values())
    result_sj = sum(len(v) for v in result_dict.values())
    compare_data = []

    # Initialize scoring system
    score_bins = np.arange(0, 1.1, 0.1)  # Score range: 0 to 1, in 0.1 increments
    score_counts = {bin: 0 for bin in score_bins}  # Total count for each bin S
    score_matches = {bin: 0 for bin in score_bins}  # Match count for each bin N

    match_sj = 0
    correct_sj = 0

    for chr_only, ref_exons in reference_dict.items():
        if chr_only in result_dict:
            res_exons = result_dict[chr_only]

            # Count total S for each score
            for res_start, res_end, score in res_exons:
                if score == 1.0:
                    score_bin = 1.0
                else:
                    score_bin = max(bin for bin in score_bins if bin <= score)
                score_counts[score_bin] += 1

            # Match reference exons
            for ref_start, ref_end in ref_exons:
                matched = False
                for res_start, res_end, score in res_exons:
                    if abs(ref_start - res_start) <= 8 and abs(ref_end - res_end) <= 8:
                        if [res_start, res_end, score] not in compare_data:
                            correct_sj += 1
                            compare_data.append([res_start, res_end, score])

                            # Update score bin matches
                            score_bin = max(bin for bin in score_bins if bin <= score)
                            score_matches[score_bin] += 1

                        matched = True
                if matched:
                    match_sj += 1



    # Calculate the match rate for each segment.
    match_percentages = {bin: (score_matches[bin] / score_counts[bin] * 100) if score_counts[bin] > 0 else 0
                         for bin in score_bins}

    # Plotting of line graphs
    plt.figure(figsize=(10, 6))
    plt.plot(score_bins, [match_percentages[bin] for bin in score_bins], marker='o', label='Match Percentage')
    plt.xlabel('Score Range')
    plt.ylabel('Match Percentage (%)')
    plt.title('Match Percentage by Score Range')
    plt.grid(True)
    plt.xticks(score_bins)
    plt.legend()
    plt.savefig('match_percentage_plot.png')
    plt.show()

    # Outputting statistics to the console
    print("Score Range | Match Percentage (%) | Total SJs | Matched SJs")
    for bin in score_bins:
        print(f"{bin:.1f}-{bin + 0.1:.1f}  | {match_percentages[bin]:>8.2f}              | {score_counts[bin]:>9} | {score_matches[bin]:>11}")

    sens = round(match_sj / reference_sj, 2) if reference_sj > 0 else 0
    prec = round(correct_sj / result_sj, 2) if correct_sj > 0 else 0

    compare_df = pd.DataFrame(compare_data, columns=['start_diff', 'end_diff', 'score'])
    compare_df.to_csv('compare.csv', index=False)

    print("sensitivity precision match_sj correct_sj reference_sj result_sj")
    print(sens, prec, match_sj, correct_sj, reference_sj, result_sj)






def main():
    args = parse_arguments()
    try:
        reference, result = load_data(args.reference, args.result)
        compare(reference, result)
    except Exception as e:
        print(f"Error occurred: {e}")


if __name__ == "__main__":
    main()

