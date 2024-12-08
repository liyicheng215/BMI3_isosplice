# l1, l2 = 
# l1: number of the vertex i (intron), the start position of i, the end position of i
# l2: the starting vertex number, the end vertex number, weight (the number of reads supporting this connection)

# import os
# print(os.getcwd())

import pandas as pd
import ast

dic = {}
SJ_pre_number = 0
gene_start_point = []
total_vertex_number = 0

def check(x, y):
    if y[7:].isdigit():
        return True
    return False

def get_vertex_number(x): # assign the vertex number
    global total_vertex_number
    if x[7:].isdigit():
        return int(x[7:])
    else:
        if x in dic.keys():
            return dic[x]
        else: # if this vertex is a begin of a gene
            total_vertex_number += 1
            dic[x] = total_vertex_number
            if x[0] == 'B': # GENE_BEGIN
                gene_start_point.append(dic[x])
            return dic[x]

def get_digit(x):
    if x[7:].isdigit():
        return int(x[7:])
    return -1

def build_graph(df):
    global SJ_pre_number
    global total_vertex_number
    for row in df.itertuples(): # get the total SJ number
        if check(row.Source, row.Target):
            # print(row.Source, row.Target)
            # print(get_digit(row.Source), get_digit(row.Target))
            SJ_pre_number = max(SJ_pre_number, max(get_digit(row.Source), get_digit(row.Target)))
    # print("SJ_pre_number = ", SJ_pre_number)
    total_vertex_number = SJ_pre_number
    graph_list = []
    for row in df.itertuples():
        if not check(row.Source, row.Target):
            continue
        source = get_vertex_number(row.Source)
        target = get_vertex_number(row.Target)
        weight = row.Weight
        if weight == '{}':
            weight_real = 0
        else:
            try: # get the weight i.e. transcript support number
                weight_dict = ast.literal_eval(weight)
                weight_real = weight_dict.get('weight', 0)
            except (ValueError, SyntaxError):
                weight_real = 0
        graph_list.append((source, target, weight_real))
    return graph_list



def graph_simplification(l2, thres, thres_dominant):
    global total_vertex_number
    n = total_vertex_number # the number of vertexs
    m = len(l2)
    print("n, m:", n, m)
    N = n + 5
    M = m + 5
    # vis = [0] * M
    ideg = [0] * N
    odeg = [0] * N
    odegmaxto = [-1] * N
    odegmax = [-1] * N
    tag = [1] * M
    vs = [[] for _ in range(N)] # vs[i] edges starting from i
    ve = [[] for _ in range(N)] # ve[i] edges ending at i
    for i in range(len(l2)): # check whether the index of the edges is starting from 0
        x, y, z = l2[i]
        # if i <= 5:
        #     print("x, y, z", x, y, z)
        # if y >= N:
        #     print("false", y)
        vs[x].append(i)
        ve[y].append(i)
        ideg[y] += z
        odeg[x] += z
        if odegmaxto[x] == -1 or z > odegmax[x]: # mark the edges starting from x with the largest weight
            odegmax[x] = z
            odegmaxto[x] = i

    # isolate the tips and bulges
    for x in range(n):
        if odeg[x] == 0 or ideg[x] == 0:
            continue
        for i in vs[x]:
            if l2[i][2] < ideg[x] * thres:
                tag[i] = 0 # tag the tips and bulges
    
    print("tagged edges = ", sum(tag[:m]))
    
    vertex_visited = [0] * N
    q = []
    for x in range(n): # check whether the index of the vertexes is starting from 0
        if ideg[x] == 0:
            q.append(x)
            vertex_visited[x] = 1
            # if not (SJ_pre_number <= x and x < total_vertex_number):
            #     print("false begin = ", x)
    
    # filtered_edges = []
    # remove all the tips and bulges
    while len(q) > 0:
        x = q.pop(0)
        for i in vs[x]:
            if tag[i] == 0:
                continue
            # filtered_edges.append(l2[i])
            # only walk through the tagged edges
            # vis[i] = 1
            y = l2[i][1]
            w = l2[i][2]
            ideg[y] -= w
            if ideg[y] == 0 and vertex_visited[y] == 0: # only go through the untagged edges
                q.append(y)
                vertex_visited[y] = 1

    # find the least confident SJs
    SJ0 = [] # the validated SJs
    for i in range(SJ_pre_number): # check whether the index of the edges is starting from 0
        if vertex_visited[i] == 0:
            SJ0.append(i)

    # find the most confident SJs
    tagSJ1 = [0] * N
    # print('N = ', N)
    for gene_start in gene_start_point:
        now = -1
        for i in vs[gene_start]:
            # print(x, now)
            x = l2[i][1]
            if now == -1 or odegmax[x] > odegmax[now]:
                now = x
        tagSJ1[now] = 1
        while 0 <= now and now < SJ_pre_number:
            # tagSJ1[now] = 1
            if odegmaxto[now] != -1 and tag[odegmaxto[now]] == 1:  # only go through the most condident edges
                if odegmax[now] >= odeg[now] * thres_dominant:
                    tag[l2[odegmaxto[now]][1]] = 1
                now = l2[odegmaxto[now]][1]
            else:
                break
    
    SJ1 = []
    for i in range(SJ_pre_number):
        if tagSJ1[i] == 1:
            SJ1.append(i)
            # if vertex_visited[i] == 0:
            #     print("false ", i)

    return SJ0, SJ1

def simplifyDAG(input_file, thres, thres_dominant):
    """Main function to simplify the Directed Acyclic Graph."""
    df = pd.read_csv(input_file)
    df.columns = ['Source', 'Target', 'Weight']

    output = build_graph(df)
    print("check SJ_pre_number = ", SJ_pre_number)
    print("total_vertex_number = ", total_vertex_number)

    SJ0, SJ1 = graph_simplification(output, thres, thres_dominant)
    print("label = 0 proportion = ", len(SJ0)/SJ_pre_number)
    print("label = 1 proportion = ", len(SJ1)/SJ_pre_number)
    print("labeled proportion = ", (len(SJ0) + len(SJ1))/SJ_pre_number)

    cluster_label = [-1] * SJ_pre_number
    # label the most confident and unconfident SJs
    for cluster in SJ0:
        cluster_label[cluster] = 0
    for cluster in SJ1:
        cluster_label[cluster] = 1
    df = pd.DataFrame({'Cluster_id': range(1, SJ_pre_number+1),
                    'Label': cluster_label})
    df.to_csv('cluster_label.csv')


if __name__ == '__main__':
    df = pd.read_csv("./edge.csv")
    df.columns = ['Source', 'Target', 'Weight']
    print(df.head())

    output = build_graph(df)
    print("check SJ_pre_number = ", SJ_pre_number)
    print("total_vertex_number = ", total_vertex_number)
    # print(len(dic))
    # print(output[0])
    # print(output[2])

    SJ0, SJ1 = graph_simplification(output, 0.01, 0.99)
    print("label = 0 proportion = ", len(SJ0)/SJ_pre_number)
    print("label = 1 proportion = ", len(SJ1)/SJ_pre_number)
    print("labeled proportion = ", (len(SJ0) + len(SJ1))/SJ_pre_number)

    cluster_label = [-1] * SJ_pre_number
    # label the most confident and unconfident SJs
    for cluster in SJ0:
        cluster_label[cluster] = 0
    for cluster in SJ1:
        cluster_label[cluster] = 1
    df = pd.DataFrame({'Cluster_id': range(1, SJ_pre_number+1),
                    'Label': cluster_label})
    df.to_csv('./Cluster_label.csv')
