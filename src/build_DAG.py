import csv
import datetime
import networkx as nx
import math
from sortedcontainers import SortedList
from collections import Counter
import numpy as np
from collections import defaultdict
 


# read in the intron information from the "find_SJ" step
def csv_read_in(file_path):
    data = []
    with open(file_path, mode='r', newline='', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.append(row)
    return data

# cluster the introns which are close to each other into one cluster and choose the one with the biggest count number as the representative
def merge_into_clusters_optimized(data, cluster_shift):
    # sort the data according to the position
    data_sorted = sorted((line[3], line[4], line[6], line[7], line[0],) for line in data)
    modified_data = data.copy() # set a new list to store the uplated data information
    cluster_info = dict()
    cluster_id = 1
    used_indices = set()  # store the indices which are already assigned
    center_map = {}  # store the information of the representative of each cluster

    for i, (chr, strain, start, end, index_id) in enumerate(data_sorted):
        if i in used_indices:
            continue
        # initialize a new cluster and record the most frequently occurring start and end positions
        cluster = {
            'index_id': [index_id],
            'starts_ends': SortedList([(start, end)]),
            'shift': 0,
            'used_indices': {i},
            'chr': chr,
            'strain': strain
        }
        # calculate and record the representative of the current cluster
        center_map[id(cluster)] = (start, end)
        modified_data[int(index_id.split("index")[1]) - 1] += ["cluster"+str(cluster_id)]  # 0-based => -1
        modified_data[int(index_id.split("index")[1]) - 1] += [0]      # shift = 0
        # merge the subsequent items into the current cluster
        j = i + 1
        while j < len(data_sorted):
            shift_score = abs(int(data_sorted[j][2]) - int(center_map[id(cluster)][0])) + abs(int(data_sorted[j][3])-int(center_map[id(cluster)][1]))
            if shift_score > cluster_shift: # compare if the shift_score smaller than our threshold
                break
            else:
                _, _, next_start, next_end, next_index_id = data_sorted[j]
                cluster['index_id'].append(next_index_id)
                cluster['starts_ends'].add((next_start, next_end))
                cluster['used_indices'].add(j)
                used_indices.add(j)
                modified_data[int(next_index_id.split("index")[1])-1] += ["cluster"+str(cluster_id)]  # 0-based => -1
                modified_data[int(next_index_id.split("index")[1])-1] += [shift_score]
                #update the representative
                counter = Counter(cluster['starts_ends'])
                most_common = counter.most_common(1)[0][0]
                center_map[id(cluster)] = most_common
                j += 1
        cluster_info["cluster"+str(cluster_id)] =[cluster['chr'], cluster['strain'], center_map[id(cluster)][0], center_map[id(cluster)][1]]
        cluster_id += 1
    return modified_data, cluster_info

# we assume that one read includes one entire gene
# cluster the reads which are close to each other into one gene and choose the one with the biggest count number as the representative
def find_gene(data, gene_shift):
    # sort the data according to the position
    data_sorted = sorted((line[3], line[4], line[5], line[1], line[0],) for line in data)
    gene_info = dict()
    gene_id = 1
    used_reads = set() # store the reads which are already assigned
    center_map = {} # store the information of the representative of each cluster
    modified_data = data.copy() 
    for i, (chr, strain, start, read_id, index_id) in enumerate(data_sorted):
        if i in used_reads:
            continue
        # initialize a new gene and record the most frequently occurring start positions
        gene = {
            'index_id': [index_id],
            'starts': SortedList([start]),
            'used_indices': {i},
            'chr': chr,
            'strain': strain
        }
        # calculate and record the representative of the current gene
        center_map[id(gene)] = start
        modified_data[int(index_id.split("index")[1])-1] += ["gene"+str(gene_id)]
        # merge the subsequent items into the current gene
        j = i + 1
        all_read_ids = [read_id]
        all_starts = [start]
        while j < len(data_sorted):
            if abs(int(data_sorted[j][2]) - int(center_map[id(gene)])) > gene_shift:
                break
            else:
                _, _, next_start, next_read_id, next_index_id = data_sorted[j]
                gene['index_id'].append(next_index_id)
                gene['starts'].add(next_start)
                gene['used_indices'].add(j)
                used_reads.add(j)
                modified_data[int(next_index_id.split("index")[1])-1] += ["gene"+str(gene_id)]
                #update the representative
                if next_read_id not in all_read_ids:
                    all_starts.append(next_start)
                counter = Counter(all_starts)
                most_common = counter.most_common(1)[0][0]
                center_map[id(gene)] = most_common
                j += 1
        gene_info["gene"+str(gene_id)] = [gene['chr'], gene['strain'],center_map[id(gene)], len(gene['index_id'])]
        gene_id += 1
    return modified_data, gene_info


def significance_of_anchor(t, genome_length):
    score = 0
    # D the maximum intron size for single anchor search
    expected_ocurrence = 1 + genome_length / (4**t)
    probability_of_fake = 1 - 1/expected_ocurrence
    eps = 10 ** (-15)
    if t >= 20:
        score += t
    else:
        score += -math.log2(eps + probability_of_fake)
    return score


def significance_of_splice(cluster_id, left_anchor, right_anchor, genome_length): # left exon length
    global splice_significance
    if left_anchor == -1: # if this is the first intron in read
        splice_significance[cluster_id] = max(splice_significance[cluster_id], significance_of_anchor(right_anchor, genome_length))
    elif right_anchor == -1: #if this is the last intron in read
        splice_significance[cluster_id] = max(splice_significance[cluster_id], significance_of_anchor(left_anchor, genome_length))
    else:
        splice_significance[cluster_id] = max(splice_significance[cluster_id], 
                                             min(significance_of_anchor(left_anchor, genome_length), significance_of_anchor(left_anchor, genome_length)))


# build the intron DAG
def build_intron_graph(intron_data, gene_info, cluster_position, genome_length):
    DAGs = []
    # for each_gene in all_gene:
    for each_gene in list(gene_info.keys()):
        row_of_gene = []
        for row in intron_data:
            if row[-1] == each_gene:
                one_row = [row[10], row[1]]
                row_of_gene.append(one_row)
        # each gene will have one DAG 
        gene_DAG = build_DAG_for_one_gene(row_of_gene, each_gene, cluster_position, genome_length)
        the_list = [(u, v, d) for u, v, d in gene_DAG.edges(data=True)]
        DAGs.extend(the_list)
    return DAGs


def build_DAG_for_one_gene(sub_data, each_gene, cluster_position, genome_length):
    global significance_of_splice
    DAG = nx.DiGraph()
    introns = dict()

    #find introns in each read
    for cluster_id, read_id in sub_data:
        _, _, begin, end = cluster_position[cluster_id]
        if read_id in introns.keys():
            introns[read_id].append([(begin, end), cluster_id])
        else:
            introns[read_id] = [[(begin, end), cluster_id]]

    for each in list(introns.keys()):
        intron = sorted(introns[each])    
        for i in range(len(intron) - 1):
            from_intron= intron[i][-1]  
            to_intron= intron[i+1][-1]
            # calculate the significance_of_splice
            if i == 0:
                significance_of_splice(int(from_intron.split("cluster")[1]), -1, abs(int(intron[i+1][0][0])-int(intron[i][0][1])), genome_length)
            if i == len(intron) - 2:
                significance_of_splice(int(to_intron.split("cluster")[1]), abs(int(intron[i+1][0][0])-int(intron[i][0][1])), -1, genome_length)
            elif i > 0:
                significance_of_splice(int(from_intron.split("cluster")[1]), abs(int(intron[i][0][0])-int(intron[i-1][0][1])), abs(int(intron[i+1][0][0])-int(intron[i][0][1])), genome_length)
            # add the edge, and store the coverage
            if not DAG.has_edge(from_intron, to_intron):
                DAG.add_edge(from_intron, to_intron, weight = 1)
            else:
                weight = DAG.edges[from_intron, to_intron]['weight'] + 1 # count the coverage of each edge
                nx.set_edge_attributes(DAG,{(from_intron, to_intron): weight}, 'weight')
        if len(intron) == 1:
            DAG.add_node(intron[0][-1])
        sources = [node for node in DAG.nodes() if DAG.in_degree(node) == 0]
        ends = [node for node in DAG.nodes() if DAG.out_degree(node) == 0]
        # add the begin and the end node for each gene
        new_source = f"BEGIN_{each_gene}"
        DAG.add_node(new_source)
        new_end = f"END_{each_gene}"
        DAG.add_node(new_end)
        for source in sources:
            DAG.add_edge(new_source, source)
        for end in ends:
            DAG.add_edge(end, new_end)
    return DAG

# calculate onlyN, mismatch, read coverage, shift score
def sj_score(clusterID, clusterInfo, cluster_intron, gene_info):
    cluster_data = cluster_intron.get(clusterID)   # [onlyN, mismatch, shift, start, end, gene]
    cluster_count = len(cluster_data)

    geneID = []
    mismatch, onlyN, shift = 0, 0, 0
    for i in cluster_data:
        onlyN += int(i[0])
        shift += int(i[2])
        if i[3] == clusterInfo[2] and i[4] == clusterInfo[3]:
            mismatch += int(i[1])
        geneID.append(i[5])

    gene_count = 0
    for gene in set(geneID):
        gene_count += gene_info.get(gene)[-1]

    read_coverage = round(cluster_count / gene_count, 2)
    mean_mismatch = round(mismatch / cluster_count, 2)
    mean_onlyN = round(onlyN / cluster_count, 2)
    mean_shift = round(shift / cluster_count, 2)

    return [read_coverage, mean_mismatch, mean_onlyN, mean_shift]


def cluster_score(intron_data, cluster_info, gene_info):
    # adjust intron_data data structure to generate a dict with clusterID key
    cluster_intron = defaultdict(list)
    for intron in intron_data:
        key = intron[-3]
        values = [intron[i] for i in [-5, -4, -2, 6, 7, -1]]   # [onlyN, mismatch, shift, start, end, gene]
        cluster_intron[key].append(values)
    # calculate score of each clusterID
    four_scores = []
    for clusterID, clusterInfo in cluster_info.items():
        info_score = [clusterID] + clusterInfo + sj_score(clusterID, clusterInfo, cluster_intron, gene_info)
        four_scores.append(info_score)
    return four_scores



def build_DAG(input_file_path, genome_length=100286401, cluster_shift=8 , gene_shift=8):
    intron_data = csv_read_in(input_file_path)
    intron_data, cluster_info = merge_into_clusters_optimized(intron_data, cluster_shift)
    intron_data, gene_info = find_gene(intron_data, gene_shift)
    cluster_scores = cluster_score(intron_data, cluster_info, gene_info)

    splice_number = len(cluster_scores)
    global splice_significance
    splice_significance = [0] * splice_number
    DAG = build_intron_graph(intron_data, gene_info, cluster_info, genome_length)
    cluster_scores = [a + [b] for a, b in zip(cluster_scores, splice_significance)]

    # delivery the result
    #np.savetxt('intron_data.csv', intron_data, delimiter=",", fmt='%s')
    column_names = ['clusterID', 'chr', "chain", 'start', 'end', 'read_coverage', 'mismatch', 'onlyN', 'shift', 'anchor']
    np.savetxt('cluster_scores.csv', cluster_scores, delimiter=",", fmt='%s', header=",".join(column_names), comments="")
    np.savetxt('DAG_edge.csv', DAG, delimiter=",", fmt='%s')

if __name__ == '__main__':
    input_file_path = "../align2intron/C_score_Introns.csv"
    print(datetime.datetime.now())
    build_DAG(input_file_path)
    print(datetime.datetime.now())