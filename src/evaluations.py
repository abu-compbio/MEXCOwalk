from utils import *
import glob
import os
import argparse
import pandas as pd
import math
import operator
from tqdm import trange, tqdm

#calculate evaluation scores
def calculate_scores_for_subnetworks(filename, edge_list, outfile,outfile_tab, num_genes):

    fhout = open(outfile, 'w')
    num_samples = len(load_patients_to_indices())  # number of patients
    gene_to_id = load_gene_to_id()
    id_to_gene = load_id_to_gene()
    data = load_gene_vs_patient_data()

    scores = []
    with open(filename) as f:
        connected_components = f.readlines()
        num_components = float(len(connected_components))
        t_coverage = 0.0
        t_mutex = 0.0
        t_density = 0.0
        t_covmutex = 0.0

        total_genes = 0
        print >>fhout, 'num_comp\tcoverage\tmutex\tcov_mutex\tdensity'
        for component_line in connected_components:
            component = component_line.split()
            coverage_list = []
            union_set = set([])
            individual_count = 0
            mult_count = 1
            for i in component:
                total_genes+=1
                if i not in data:
                    print(i + " not found" + filename)
                    continue
                union_set = set(data[i]).union(union_set)
                individual_count = individual_count + len(data[i])
                mult_count *= (float(len(data[i]))/float(num_samples))
                coverage_list.append(len(data[i]) / float(num_samples))
            union_set = len(list(union_set))#?

            num_in = 0
            num_out = 0
            for edge in edge_list:
                a = id_to_gene[edge[0]] in component
                b = id_to_gene[edge[1]] in component
                if a and b:
                    num_in = num_in + 1
                elif a or b:
                    num_out = num_out + 1

            density = 0.0
            if num_in + num_out != 0:
                density = round(float(num_in) / (num_in + num_out),3)

            try:
                coverage = round(float(union_set) / num_samples,3)
                mutex = round(float(union_set) / individual_count,3) if individual_count != 0 else 0
            except:
                print('file: ',outfile)
            cov_mutex = mutex*coverage

            t_density = t_density + density
            t_coverage = t_coverage + coverage
            t_mutex = t_mutex + mutex
            t_covmutex = t_covmutex + cov_mutex

            print >>fhout, str(len(component)) + '\t' + str(coverage) +\
            '\t' + str(mutex) + '\t' +  str(cov_mutex) + '\t' +\
            str(density)

            scores.append([len(component), coverage, mutex, cov_mutex, density])

    avg_density = t_density / num_components
    avg_cov = t_coverage / num_components
    avg_mutex = t_mutex / num_components
    avg_covmutex = t_covmutex / num_components


    print >>fhout, '****AVERAGE ACROSS ALL MODULES*****'
    print >>fhout, str(total_genes) + '\t' + str(int(num_components)) + \
    '\t' + str(t_density) + '\t' + str(avg_density) + \
    '\t' + str(t_coverage) + '\t' + str(avg_cov) + \
    '\t' + str(t_mutex) + '\t' + str(avg_mutex) + \
    '\t' + str(t_covmutex) + '\t' + str(avg_covmutex)
    fhout.close()

    fhout = open(outfile_tab, 'a+')
    print >>fhout, str(num_genes) + '\t' + str(int(num_components)) + \
    '\t' + str(t_density) + '\t' + str(avg_density) + \
    '\t' + str(t_coverage) + '\t' + str(avg_cov) + \
    '\t' + str(t_mutex) + '\t' + str(avg_mutex) + \
    '\t' + str(t_covmutex) + '\t' + str(avg_covmutex)
    return [t_density, avg_density, t_coverage, avg_cov, t_mutex, avg_mutex, t_covmutex, avg_covmutex]   #[t_coverage, t_mutex, t_density, t_density+t_coverage+t_mutex]

#read genes
def read_all_genes(filename):
    l = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            for i in line:
                l.append(i)
    return l

#****** COSMIC RELATED*******
def read_subnetworks(filename):
    fhmodule = open(filename)
    module_genes = []
    module_list = []
    huge_module = []
    for line in fhmodule:
       words = line.strip().split()
       module_list.append(words)
       for word in words:
           if word not in module_genes:
              module_genes.append(word)
       if 'TP53' in words:
           huge_module = words

    return (module_genes, module_list, huge_module)

def find_overlap(module, cosmic_genes):
    count = 0
    for gene in module:
       if gene in cosmic_genes:
           count += 1
    return (count, len(module))

#find number of genes that overlap with cosmic genes
def cosmic_overlap_analysis(file, cosmic_genes, module_name, num_genes):
    #print file

    cosmic_file = '../out/evaluation/{}.txt'.format(module_name)
    cosmic_huge_file = '../out/evaluation/{}_huge.txt'.format(module_name)
    cosmic_tab_file = '../out/evaluation_tab/{}.txt'.format(module_name)

    largest_module = ''
    largest_module_cosmic = ''
    x_axes = ''
    (module_genes, module_list, huge_module) = read_subnetworks(file)
    (count_overlap, count) =  find_overlap(module_genes, cosmic_genes)

    with open(cosmic_file, "a") as f:
        f.write("***Top " + str(num_genes) + " genes***\n")
        f.write("Overlap with cosmic: " + str(count_overlap) + "\n")
        count_overlap_huge = find_overlap(huge_module, cosmic_genes)[0]
        f.write("Overlap with cosmic only huge module: "  + str(count_overlap_huge) +  "\n\n" )
        x_axes += str(num_genes) + '\t'
        largest_module += str(len(huge_module)) + '\t'
        largest_module_cosmic += str(count_overlap_huge) + '\t'

    with open(cosmic_huge_file, "w+") as f:
      f.write(x_axes + '\n')
      f.write(largest_module + '\n')
      f.write(largest_module_cosmic + '\n')


    with open(cosmic_tab_file, "a") as f:
      # number of genes, overlap
      f.write(str(num_genes) + "\t" + str(count_overlap) + "\n")

    return count_overlap

#generate cosmic overlap file for evaluation
def generate_cosmic_analysis_file(path_pre = "../out/evaluation_tab/"):
    d = {}
    for key in tqdm(models, desc='running cosmic analysis'):


        with open(path_pre + key + "/cosmic_" + key + ".txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                N = line[0].strip()
                our = line[1].strip()
                d[str(N) + key] = our

    with open(path_pre + "cosmic" + ".txt", "w") as ff:
        ff.write('N\t')
        for k in models:
            ff.write(k + '\t')
        ff.write('\n')
        num_genes_list = range(ne,ns+step, step)
        for i in num_genes_list:
            ff.write(str(i))
            for key in models:
                if str(i) + key in d:
                    ff.write('\t' + d[str(i) + key])
                else:
                    ff.write('\t' + ' ')
            ff.write('\n')


#generate all evaluation files
def generate_eval_files(path_pre = "../out/evaluation_tab/"):
    eval_list = ['num_components', \
    't_density', 'avg_density', \
    't_coverage', 'avg_cov', \
    't_mutex', 'avg_mutex', \
    't_covmutex', 'avg_covmutex']

    for i in trange(len(eval_list), desc='running evaluation'):
        eval = eval_list[i]
        d = {}
        for key in models:

            with open(path_pre + key + "/eval_" + key + ".txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.split()
                    N = line[0].strip()
                    e = line[i+1].strip()
                    d[str(N) + key] = e


        with open(path_pre + eval + ".txt", "w") as f:
            f.write('N\t')
            for k in models:
                f.write(k + '\t')
            f.write('\n')
            num_genes_list = range(ne,ns+step, step)
            for i in num_genes_list:
                f.write(str(i))
                for key in models:
                    if str(i) + key in d:
                        f.write('\t' + d[str(i) + key])
                    else:
                        f.write('\t' + ' ')
                f.write('\n')

# calculate weighted mutex, coverage and density scores
def calculate_weighted_scores():
    evals = [ 'iwavg_cov', 'wavg_mutex', 'wavg_density']
    evals_pos = [1, 2, 4]

    dir_pre = "../out/evaluation/"
    dir_post = "/optimized_function_comparison/"


    for i in trange(len(evals), desc='Calculating weighte scores'):
        eval = evals[i]
        pos = evals_pos[i]
        with open("../out/evaluation_tab/"+ eval + ".txt", "w") as eval_file:
            eval_file.write("N\t")
            for model in models:
                eval_file.write(model + "\t")

            eval_file.write("\n")
            num_genes_list =range(ne,ns+step, step)

            for n in num_genes_list:
                eval_file.write(str(n) + "\t")
                for model in models:
                    score = 0.0
                    filepath = dir_pre + model + dir_post + model + "_" + str(n) + ".txt"

                    if os.path.exists(filepath):
                        with open(dir_pre + model + dir_post + model + "_" + str(n) + ".txt") as n_file:
                            lines = n_file.readlines()[1:-2]
                            if pos == 1:
                                score = 0.0
                                genes = 0
                                weight_sum = 0
                                for line in lines:
                                    line = line.strip().split()
                                    weight_sum += (1-(float(line[0])/n))
                                    score += ((1-(float(line[0])/n))) * float(line[pos])

                                    genes += float(line[0])
                                score /= weight_sum
                            else:
                                score = 0.0
                                genes = 0
                                for line in lines:
                                    line = line.strip().split()
                                    score += float(line[0])*float(line[pos])
                                    genes += float(line[0])
                                score /= genes
                        eval_file.write(str(score) + '\t')
                    else:
                        eval_file.write(' ' + '\t')

                eval_file.write("\n")

#prepare file paths into a dictionary with model as key
def prep_file_paths(key):
    paths_raw = glob.glob(key)

    dict_paths = {}
    for filename in paths_raw:

          directory, file_core = os.path.split(filename)

          if file_core.startswith('cc'):
              num_genes = int(file_core.split('_')[1][1:])
              dict_paths[num_genes] = filename

    sorted_dict = sorted(dict_paths.items(), key=operator.itemgetter(0))

    paths = []
    for i in range(len(sorted_dict)):
        paths.append(sorted_dict[i][1])
    return paths

#calculate inverse average coverage*average mutex
def calculate_product_iwcov_wmex():
    filepath_iwcov = "../out/evaluation_tab/iwavg_cov.txt"
    filepath_wmex = "../out/evaluation_tab/wavg_mutex.txt"
    filepath_iwcov_wmex = "../out/evaluation_tab/iwcov_wmex.txt"

    top_genes_to_read = args.numstart - args.numend + args.stepsize

    iwcov_df = pd.read_table(filepath_iwcov, nrows=top_genes_to_read / args.stepsize, index_col=0, na_values=' ')
    wmex_df = pd.read_table(filepath_wmex, nrows=top_genes_to_read / args.stepsize, index_col=0, na_values=' ')

    iwcov_df.drop(iwcov_df.columns[len(iwcov_df.columns) - 1], axis=1, inplace=True)
    wmex_df.drop(wmex_df.columns[len(wmex_df.columns) - 1], axis=1, inplace=True)

    mult = iwcov_df * wmex_df
    mult.to_csv(filepath_iwcov_wmex, sep='\t')


if __name__ == "__main__":

    # parse arguments
    description = "Evaluate the model(s) and calculate scores"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', "--models", type=str, required=False, nargs='+',  default=['mexcowalk'],
                        help="model(s) to evaluate")
    parser.add_argument('-ns', '--numstart', type=int, required=False, default=1000, help="max number of genes")
    parser.add_argument('-ne', '--numend', type=int, required=False, default=100, help="min number of genes")
    parser.add_argument('-step', '--stepsize', type=int, required=False, default=100, \
                        help="step size of decrement from num start to num end")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose')

    args = parser.parse_args()


    models = args.models
    ns = args.numstart
    ne = args.numend
    step = args.stepsize
    print models
    # Read the file once, instead of reading it in a loop many times
    cosmic_genes = read_cosmic_genes()

    #evaluate all models provided
    for key in tqdm(models):

        current_path_key = key

        paths = prep_file_paths('../out/connected_components_isolarge/' + current_path_key +'/*.txt')


        newpath = '../out/evaluation/{0}'.format(current_path_key)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            newpathext = newpath + '/optimized_function_comparison'
            if not os.path.exists(newpathext):
                os.makedirs(newpathext)

        newpath = '../out/evaluation_tab/{0}'.format(current_path_key)
        if not os.path.exists(newpath):
            os.makedirs(newpath)


        our_list = []
        edge_list_original = load_edge_list()
        fhout = open('../out/evaluation/'+ current_path_key + '/optimized_function_comparison/summary_' + current_path_key  +'.txt', 'w')
        num_genes_list = range(ne,ns+step, step)

        #initializing cosmic overlap files for the specific model
        cosm_file = open('../out/evaluation/{}/cosmic_{}.txt'.format(key,key), 'w')
        cosm_tab_file = open('../out/evaluation_tab/{}/cosmic_{}.txt'.format(key, key), 'w')
        eval_tab_file = open('../out/evaluation_tab/{}/eval_{}.txt'.format(key, key), 'w')
        cosm_file.close()
        cosm_tab_file.close()
        eval_tab_file.close()

        for i in trange(min(((ns-ne+step)/step), len(paths)), desc='running main evaluation'):

            num_genes = str(num_genes_list[i])

            current_path = glob.glob('../out/connected_components_isolarge/{}/cc_n{}_*'.format(key,num_genes))[0]

            module_name = current_path_key + '/cosmic_' + current_path_key
            overlap = cosmic_overlap_analysis(current_path, cosmic_genes, module_name, num_genes)
            print >>fhout, overlap,

            outfile = '../out/evaluation/'+ current_path_key + '/optimized_function_comparison/' + current_path_key + '_' + str(num_genes) + '.txt'
            outfile_tab = '../out/evaluation_tab/'+ current_path_key + '/eval_' + current_path_key + '.txt'

            scores = calculate_scores_for_subnetworks(current_path, edge_list_original, outfile, outfile_tab, num_genes)

            for j in range(len(scores)):
                print >>fhout, scores[j],
            print >>fhout, '\n',


        fhout.close()

    generate_cosmic_analysis_file()
    generate_eval_files()
    calculate_weighted_scores()
    calculate_product_iwcov_wmex()
