#!/usr/bin/env python

import networkx as nx
from utils import *
import numpy as np
import scipy.sparse as sp
import operator
import os
from collections import deque
import glob
import argparse

#find the largest component from a given list of graphs
def find_largest_comp(list_comp_graph):
    large = list_comp_graph[0]
    large_index = 0
    for i in range(len(list_comp_graph)):
        if len(large) < len(list_comp_graph[i]):
            large = list_comp_graph[i]
            large_index = i

    return (large,large_index)

#convert list of graphs to list of components
def list_graph_to_list_comp(list_graph):
    list_comp = []
    for i in range (len(list_graph)):
        list_comp.extend([list_graph[i].nodes])

    return list_comp

#find largest component length
def max_comp_length(list_graph):
    max_length = len(list_graph[0])

    for i in range (len(list_graph)):
        if len(list_graph[i])> max_length:
            max_length = len(list_graph[i])

    return max_length

#find gene with max out degree wuthin the constructed graph
def find_max_outdegree(comp):
    max_out_degree = 0
    largest_node_id = 0
    largest_gene_id = ''

    for n in comp.nodes():
        out_n = comp.out_degree(n)
        if max_out_degree < out_n:
            max_out_degree = out_n
            largest_node_id = n
            largest_gene_id = id_to_gene[n]

    return (max_out_degree, largest_node_id, largest_gene_id)

#star construction by finding all interconnected genes
def star_construction(large_comp, largest_node_id):
    all_neighbor_of_largest = set()
    star_nodes_set = set([largest_node_id])

    for e in large_comp.edges:

        if largest_node_id in e:

            all_neighbor_of_largest.update(e)

    all_neighbor_of_largest.remove(largest_node_id)

    for u in all_neighbor_of_largest:
        in_comp = True
        for e in large_comp.edges:
            if e[0] == u:
                v = e[1]

                if v != largest_node_id and (v not in all_neighbor_of_largest):

                    in_comp = False

                    break
            elif e[1] == u:
                v = e[0]

                if v != largest_node_id and (v not in all_neighbor_of_largest):
                    in_comp = False

                    break

        if in_comp == True:
            star_nodes_set.update([u])

    return list(star_nodes_set)



#evaluates each component and returns a split component graph for largest components
def split_large_components(list_comp_graph, large_threshold = 0):


    list_comp = list_graph_to_list_comp(list_comp_graph)
    set_all_genes_before = set([n for comp in list_comp for n in comp])

    threshold_small_comp = args.k_threshold


    gene_set = set()
    for i in range(len(list_comp)):
        gene_set = gene_set.union(set(list_comp[i]))

    max_out_degree_all = 0

    #finding large graph threshold, take the max of max-out-degree of each component
    for comp in list_comp_graph:
         max_out_degree_pre, largest_node_id_pre, largest_gene_id_pre = find_max_outdegree(comp)

         if max_out_degree_all < max_out_degree_pre:
               max_out_degree_all = max_out_degree_pre


    if large_threshold == 0:
        threshold_large_comp = max_out_degree_all
    else:
        threshold_large_comp = large_threshold


    list_graph_leftover = []
    list_large_components = []
    for comp in list_comp_graph:
        if len(comp)>= threshold_large_comp:
            list_large_components.append(comp)

        else:
            list_graph_leftover.append(comp)

    #going through all large components
    all_modified_component_list = []

    for lc in list_large_components:
        main_comp_list = []
        small_comp_list = []
        #large_graph_queue
        large_comp_queue = deque()
        large_comp_queue.append(lc)
        while len(large_comp_queue) >0:

            largest_comp = large_comp_queue.popleft()
            max_out_degree, largest_node_id, largest_gene_id = find_max_outdegree(largest_comp)
            reduced_comps = largest_comp.copy()

            removable_nodes_list = star_construction(largest_comp,largest_node_id)

            #reduced comps -> largest graph - large gene and neighbors
            reduced_comps.remove_nodes_from(removable_nodes_list)

            #adding to LC or SCL
            #checking if star construction is less than k
            largest_gene_graph = largest_comp.subgraph(removable_nodes_list).copy()
            if len(largest_gene_graph) < threshold_small_comp:
                small_comp_list.append(largest_gene_graph)
            else:
                main_comp_list.append(largest_gene_graph)

            scc = nx.strongly_connected_components(reduced_comps)

            for comp in scc:

                if len(comp)> threshold_large_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    large_comp_queue.append(g)

                elif len(comp)< threshold_small_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    small_comp_list.append(g)

                else:
                    g = reduced_comps.subgraph(comp).copy()
                    main_comp_list.append(g)


        temp_list = []
        for scm in range(len(main_comp_list)):

            if len(main_comp_list[scm]) < threshold_small_comp:
                small_comp_list.append(main_comp_list[scm])

            else:
                temp_list.append(main_comp_list[scm])

        main_comp_list = temp_list[:]

        #adding small components less than k back into one of the components based on their scores
        for scomp in range(len(small_comp_list)):
            max_comp_score = 0
            max_comp_index = 0
            for comp_index in  range(len(main_comp_list)):
                comp_score = 0
                for gene_m in main_comp_list[comp_index]:

                    for gene_s in small_comp_list[scomp]:

                        if (gene_s, gene_m) in lc.edges:
                            comp_score +=1
                        if (gene_m, gene_s) in lc.edges:
                            comp_score += 1

                if comp_score > max_comp_score:
                    max_comp_score = comp_score
                    max_comp_index = comp_index

            main_comp_list[max_comp_index].add_nodes_from(small_comp_list[scomp])

            temp_subgraph_nodes = main_comp_list[max_comp_index].nodes

            # also add the edges
            main_comp_list[max_comp_index] = lc.subgraph(temp_subgraph_nodes).copy()


        all_modified_component_list.extend(main_comp_list[:])

    set_all_genes_after = set([n for comp in list_graph_leftover[:] + all_modified_component_list[:] for n in comp.nodes])


    assert set_all_genes_after == set_all_genes_before
    return list_graph_leftover[:] + all_modified_component_list[:]


if __name__ == '__main__':

    # parse arguments
    #set starting threshold at 0.00008 and set numstart to 2500 for top 2500 genes
    #default starting threshold at 0.0002 for top 1000 genes

    description = "Find connected components then isolate and split large components"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', "--model", type=str, required=False, default='mexcowalk',
                        help="model to find connected components")
    parser.add_argument('-k', "--k_threshold", type=int, required=False, default=3,
                        help="min length of module")
    parser.add_argument('-ns', '--numstart', type=int, required=False, default=1000, help="max number of genes")
    parser.add_argument('-ne', '--numend', type=int, required=False, default=100, help="min number of genes")
    parser.add_argument('-step', '--stepsize', type=int, required=False, default=100, \
                        help="step size of decrement from num start to num end")
    parser.add_argument('-ts', '--start_threshold', type=float, required=False, default=0.0002, help="starting threshold")
    parser.add_argument('-cc', '--cc_original', action='store_true', \
                        default=False, help="generate connected component files pre-split")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose')

    args = parser.parse_args()

    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id()

    model = args.model

    #input file path
    filepath = "../out/connected_components_original/" + model + "/"
    cc_path = "../out/connected_components_isolarge/" + model + "/"

    try:
        # Get the starting threshold from the file name at n = numstart
        threshold_start = float(glob.glob(filepath+'cc_n{}_*'.format(args.numstart))[0].split('/')[-1].split('d')[1][:-4])-0.000001
    except:
        threshold_start = args.start_threshold

    path = cc_path

    if not os.path.exists(cc_path):
        os.makedirs(cc_path)
    if not os.path.exists(filepath) and args.cc_original:
        os.makedirs(filepath)

    LARGE_NUM = 100000
    k = args.k_threshold
    our_E = sp.load_npz("../out/random_walk/"+ model+"_sparse_matrix_e.npz")
    E = our_E.toarray()
    # find threshold
    num_start = args.numstart  # max n
    num_end = args.numend  # min n

    # create the initial graph, omit the edges smaller than the threshold

    N = len(id_to_gene)
    G = nx.DiGraph()

    for i in range(N):
        G.add_node(i)

    for i in range(N):
        for j in range(N):
            if E[i][j] > threshold_start:
                G.add_edge(j, i)

    # find the list of the connected components
    list_graphs = []
    num_nodes = 0
    count_comp = 0
    smallest_weight = LARGE_NUM
    smallest_edge_str = ''
    for comp in nx.strongly_connected_components(G):
        if len(comp) >= k:
            subG = G.subgraph(comp).copy()
            list_graphs.append(subG)
            num_nodes += len(comp)

            for e in subG.edges():
                u = e[0]
                v = e[1]
                key = str(count_comp) + '_' + str(u) + '_' + str(v)
                if E[v][u] < smallest_weight:
                    smallest_edge_str = key
                    smallest_weight = E[v][u]
            count_comp += 1
            if args.verbose:
                print smallest_edge_str

    if args.verbose:
        print path
        print 'number of nodes in the beginning is : ', num_nodes
        print smallest_edge_str

    while num_start >= num_end:
        print("a" + str(num_start))
        iteration = 0
        threshold_found = 0
        num_target_genes = num_start
        while num_nodes > num_target_genes:
            if args.verbose:
                print(str(num_target_genes) + " " + str(num_nodes) + " " + str(iteration))

            # remove smallest edge

            smallest_edge = smallest_edge_str.split('_')
            threshold_found = smallest_weight
            graph_index = int(smallest_edge[0])

            subG = list_graphs[graph_index]
            num_nodes -= len(subG.nodes())


            subG.remove_edge(int(smallest_edge[1]), int(smallest_edge[2]))
            subG_comps = nx.strongly_connected_components(subG)

            del list_graphs[graph_index]
            for comp in subG_comps:
                if len(comp) >= k:
                    list_graphs.append(subG.subgraph(comp).copy())
                    num_nodes += len(comp)

            num_comps = len(list_graphs)

            smallest_weight = LARGE_NUM
            smallest_edge_str = ''
            for i in range(num_comps):
                graph = list_graphs[i]

                if len(graph.nodes()) >= k:
                    for e in graph.edges():
                        u = e[0]
                        v = e[1]
                        if E[v][u] < smallest_weight:
                            key = str(i) + '_' + str(u) + '_' + str(v)
                            smallest_edge_str = key
                            smallest_weight = E[v][u]
            iteration += 1

        if args.cc_original:
            # filename contains number of genes and the threshold for obtaining these number of genes respectively
            outfilename = filepath + "cc_n" + str(num_target_genes) + "_d" + str(threshold_found) + ".txt"


            #prints the top [largest component size - difference] genes (or target genes)
            with open(outfilename, "w") as f:
                for i in range(len(list_graphs)):
                    for node_index in list_graphs[i]:
                        f.write(id_to_gene[node_index] + " ")
                    f.write("\n")
                f.close()

        ################################################################################################################
        #isolating and splitting largest components

        largest_comp, largest_comp_index = find_largest_comp(list_graphs)
        component_list = list_graphs[:]
        largeset = set([id_to_gene[idx] for idx in largest_comp])

        largest_gene_comp_list = []

        #inputs = component list, threshold (0 means largest component out degree)
        final_comp_list_graphs = split_large_components(component_list)

        # filename contains number of genes, length of the largest component, number of connected components
        # and the threshold for obtaining these number of genes respectively
        outfilename = path + "cc_n" + str(num_target_genes) + "_"+ str(len(largest_comp))+ "_"+\
                      str(len(final_comp_list_graphs)) + "_d" + str(threshold_found) + ".txt"


        #prints the top [largest component size - difference] genes (or target genes)
        with open(outfilename, "w") as f:
            for i in range(len(final_comp_list_graphs)):
                for node_index in final_comp_list_graphs[i]:
                    f.write(id_to_gene[node_index] + " ")
                f.write("\n")
            f.close()


        num_start -= args.stepsize
