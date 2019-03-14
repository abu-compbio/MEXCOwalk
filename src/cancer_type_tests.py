import scipy.stats as stats
import math
import os
import glob
# 'ClusterOne_0_0.9_0.6_k3',
models = [
        "hotnet2",
        "memcover_v1","memcover_v2","memcover_v3",
        "mutex_t07_nsep_cov",
        'hier_hotnet'
        ]


num_subtypes = 11 #combine CRC
subtypes = ['BLCA', 'BRCA', 'CRC', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC' ]

def count_subtype(patient_list):
   count_vec = [0] *  num_subtypes
   for patient in patient_list:
      subtype = patient_types[patient]
      index = subtypes.index(subtype)
      count_vec[index] += 1

   return count_vec



fhmutation = open("../data/genes_vs_patients_indices_gen_paran.txt")

fhindices = open("../data/patients_to_indices_gen.txt")

ftype = open("../../memcover/data/our_patients_subtypes.txt")


# read the patients / types file

count_types_each = [0] * num_subtypes
patient_types = {}
for line in ftype:
	words = line.strip().split()
	patient_types[words[0]] = words[1]
	subtype_index = subtypes.index(words[1])
	count_types_each[subtype_index] += 1

print (count_types_each)
count_subtypes_except_this = [0] * num_subtypes
for i in range(num_subtypes):
	count_subtypes_except_this[i] = sum(count_types_each) - count_types_each[i]

print ('sum(count_types_each): ',sum(count_types_each))
print ('count_subtypes_except_this: ',count_subtypes_except_this)

dict_indices = {}
for line in fhindices:
	words = line.strip().split()
	dict_indices[words[1]] = words[0]

dict_genes_to_indices = {}
dict_genes_to_patients = {}
for line in fhmutation:
	words = line.strip().split()
	gene = words[0]
	indices = words[1:]
	dict_genes_to_indices[gene] = indices
	dict_genes_to_patients[gene] = []
	for index in indices:
		dict_genes_to_patients[gene].append(dict_indices[index])


#entry_11 entry_12
#entry_21 entry_22


list_range = list(range(100, 600, 100))
list_range.append(554)
list_range.extend([600,700,800])
list_range.append(806)
list_range.extend(list(range(900,1601,100)))
memcover_v3_range = list_range[:]
list_range.extend(list(range(1700,2501,100)))
print (list_range)


common_path = '../hint/out/connected_components_isolarge_n2500_whh/'
output_path = '../hint/out/cancer_subtype_test_results_2/'


fhout_signif_per_comp = open(output_path + 'Signif_tests_per_comp.txt', 'w')
fhout_logpval_per_comp = open(output_path + 'Logpval_per_comp.txt', 'w')

dict_signif = {}
dict_logpval = {}

for key in models:
    # os.system('mkdir ' + output_path + '/' + key)
    if not os.path.exists(output_path  + '/' + key):
        os.mkdir(output_path  + '/' + key)
    dict_signif[key] = {}
    dict_logpval[key] = {}
    path = common_path + key + '/'


    for N in list_range:
        skip = False
        if key == 'memcover_v3' and N > 1600:
        	skip = True
        if key == 'hier_hotnet' and N not in [554, 806]:
        	skip = True

        if skip == False:
            print ('path is ', path + 'cc_n' + str(N) + '_*')
            comp_filename = glob.glob(path+'cc_n' + str(N) + '_*')[0]
            # print 'comp filename ' ,comp_filename
            fhcomponents = open(comp_filename)
            #fhcomponents = open("../hint/out/connected_components_isolarge_v2/mutex_t07_nsep_cov_nsep/cc_n100_48_10_d0.000327361545892.txt")
            #fhcomponents = open("../../memcover/results/components/cc_n100_memcover_hint_0.548_k0_genes.txt")

            #fhcomponents = open("../../memcover/results/components/cc_n100_memcover_hint_0.03_k0_genes.txt")

            #fhcomponents = open("../hint/out/connected_components_isolarge_v2/hotnet2/cc_n100_k3_d0.000395306627809.txt")



            #fhout = open("../hint/out/cc_n100_memcover_hint_0.548_k0_genes_cancer_subtype_tests.txt", 'w')

            fhout = open(output_path + '/' + key + '/cc_n' + str(N)  + '_cancer_subtype_tests.txt'  , 'w')



            comp_index = 0

            num_components = len(fhcomponents.readlines())
            fhcomponents.seek(0)

            subtype_coverage = [0] * num_subtypes

            sum_logpval = 0
            min_logpvals = []
            for component in fhcomponents:
                genes = component.strip().split()
                #print component
                patient_list_union = []
                for gene in genes:
                	patient_list_union.extend(dict_genes_to_patients[gene])
                #print 'not unique ', len(patient_list_union)
                patient_list_union = list(set(patient_list_union))
                #print 'patient list union unique ', len(patient_list_union)
                count_mutated_subtypes_vec = count_subtype(patient_list_union)
                fhout.write('COMPONENT is {} {}\n'.format(comp_index, genes))
                # print >>fhout, 'COMPONENT is ', comp_index, genes
                min_pval = 1
                for i in range(num_subtypes):
                    entry_11 = count_mutated_subtypes_vec[i]
                    entry_12 = count_types_each[i] - count_mutated_subtypes_vec[i]
                    entry_21 = sum(count_mutated_subtypes_vec) - count_mutated_subtypes_vec[i]
                    entry_22 = count_subtypes_except_this[i] - entry_21
                    oddsratio, pvalue = stats.fisher_exact([[entry_11, entry_12], [entry_21, entry_22]])
                    # print >>fhout, 'subtype ', subtypes[i] , entry_11, entry_12, entry_21, entry_22, pvalue,
                    corrected_pvalue = pvalue * num_components * num_subtypes
                    logpval = math.log(corrected_pvalue, 10)
                    fhout.write('subtype is {} {}'.format(corrected_pvalue, round(logpval,2)))
                    # print >> fhout, 'subtype', corrected_pvalue, round(logpval,2),
                    if logpval <= min_pval:
                    	min_pval = logpval

                    if corrected_pvalue < 0.05:
                        fhout.write(' ***SIGNIFICANT***\n')
                        # print >>fhout, '***SIGNIFICANT***'
                    	subtype_coverage[i] += 1
                    else:
                    	# print >>fhout, ''
                        fhout.write('\n')

                min_logpvals.append(min_pval)
                sum_logpval += min_pval
                comp_index += 1

            fhout.write(' '.join([str(xyz) for xyz in  min_logpvals]))
            print (subtype_coverage)
            print (sum(subtype_coverage) / float(num_components))
            print (sum_logpval / float(num_components))
            dict_signif[key][N] = round(sum(subtype_coverage) / float(num_components),2)
            dict_logpval[key][N] = round(sum_logpval / float(num_components),2)
            fhout.close()

fhout_signif_per_comp.write('\t'.join([str(elem) for elem in list_range]))
fhout_logpval_per_comp.write('\t'.join([str(elem) for elem in list_range]))
# print >>fhout_signif_per_comp, '\t'.join([str(elem) for elem in list_range])
# print >>fhout_logpval_per_comp, '\t'.join([str(elem) for elem in list_range])
for key in models:
    forprint_signif = key
    forprint_logpval = key
    for N in list_range:
        if (N not in [554, 806] and key == 'hier_hotnet') or (N not in memcover_v3_range and key == 'memcover_v3'):
            forprint_signif += '\tN/A'
            forprint_logpval += '\tN/A'
        else:
            forprint_signif += '\t' + str(dict_signif[key][N])
            forprint_logpval += '\t' + str(dict_logpval[key][N])

    fhout_signif_per_comp.write(forprint_signif+'\n')
    fhout_logpval_per_comp.write(forprint_logpval+'\n')

    # print >>fhout_signif_per_comp, forprint_signif
	# print >>fhout_logpval_per_comp,	forprint_logpval
