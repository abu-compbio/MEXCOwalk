#!/usr/bin/env python

from utils import *
import scipy.sparse as sp
import numpy as np
import os
import argparse

def normalizeW(A):
    n = len(A)
    W = np.zeros((n, n))
    for j in range(n):
        d_j = float(A[:,j].sum())
        if(d_j == 0):    # added this check
            continue
        for i in range(n):
            W[i,j] = A[i,j]/d_j
    return W

def computeF(W, beta):
    n = len(W)
    return beta*np.linalg.inv(np.eye(n)-(1-beta)*W)

def computeH_pan():
    gene_to_id = load_gene_to_id()
    lines = []
    with open("../data/" + gene_score_file) as f:
        lines = f.readlines()
    N = len(gene_to_id)
    h = np.zeros((N, N))
    yes = 0
    no = 0
    for line in lines:
        line = line.split()
        if line[0] not in gene_to_id.keys():
            no += 1
            continue
        yes += 1
        gene_id = gene_to_id[line[0]]
        h[gene_id][gene_id] = line[1]
    #print(yes)
    #print(no)
    return h

def computeH():
    # some genes in patient_data are not in the graph, so we have to ignore those
    data = load_gene_vs_patient_data()
    N = len(load_patients_to_indices())  # number of patients
    gene_list = load_gene_list()
    genes = {}

    for key in gene_list:
        if key in data:
            genes[key] = len(data[key]) / float(N)

    N = len(gene_list) # number of genes
    h = np.zeros((N, N))
    gene_to_id = load_gene_to_id()
    no = 0
    for key in genes:
        # TODO there are genes in patiend data that does not exist in the graph
        if key not in gene_to_id:
            no += 1
            continue
        gene_id = gene_to_id[key]
        h[gene_id][gene_id] = genes[key]
    #print(no)
    return h

def computeW():
    gene_to_id = load_gene_to_id()

    N = len(gene_to_id)
    w = np.zeros((N, N))

    with open(weight_out_file) as f:
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = line.split()
            w[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
    return normalizeW(w)

def computeW_combine(alpha1, alpha2):
    gene_to_id = load_gene_to_id()

    N = len(gene_to_id)
    w1 = np.zeros((N, N))
    w2 = np.zeros((N, N))
    with open(weight_out_file) as f:
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = line.split()
            #print line
            w1[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w1[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
            w2[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[3]
            w2[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[3]


    w1_norm = normalizeW(w1)
    w2_norm = normalizeW(w1)
    return alpha1 * w1_norm + alpha2 * w2_norm

if __name__ == "__main__":
    network_file = 'hint'

    # parse arguments
    description = "Perform random walk on edge weighted graph"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', "--model", type=str, required=False, default='mexcowalk', help="model to perform random walk")
    parser.add_argument('-s', '--store', action='store_true', default=False, help="store all matrices")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose')

    args = parser.parse_args()
    model = args.model

    path_pre = "../out/random_walk/"

    if not os.path.exists(path_pre):
        os.makedirs(path_pre)

    weight_out_file = "../out/edge_weights/" + model + '.txt'

    W = computeW()
    if args.store:
        sp_w = sp.csc_matrix(W)
        w_path = path_pre + model + "_sparse_matrix_w.npz"
        sp.save_npz(w_path, sp_w)
    if args.verbose:
        print("W")
        print(W)

    F = computeF(W, network_beta)
    if args.store:
        sp_f = sp.csc_matrix(F)
        f_path = path_pre + model + "_sparse_matrix_f.npz"
        sp.save_npz(f_path, sp_f)
    if args.verbose:
        print("F")
        print(F)

    H = computeH_pan()
    if args.store:
        sp_h = sp.csc_matrix(H)
        h_path = path_pre  + model + "_sparse_matrix_h.npz"
        sp.save_npz(h_path, sp_h)
    if args.verbose:
        print("H")
        print(H)

    E = np.dot(F, H)
    sp_e = sp.csc_matrix(E)
    e_path = path_pre  + model +"_sparse_matrix_e.npz"
    sp.save_npz(e_path, sp_e)
    if args.verbose:
        print("E")
        print(E)

