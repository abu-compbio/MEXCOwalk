import scipy.stats as stats

num_subtypes = 11 #combine CRC

subtypes = ['BLCA', 'BRCA', 'CRC', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC' ]

def count_subtype(patient_list):
   count_vec = [0] *  num_subtypes
   for patient in patient_list:
      subtype = patient_types[patient]
      index = subtypes.index(subtype)
      count_vec[index] += 1
      
   return count_vec 

   


fhcomponents = open()

fhmutation = open("../data/genes_vs_patients_indices_gen_paran.txt")

fhindices = open("../data/patients_to_indices_gen.txt")

fhout = open("")

# read the patients / types file

count_types_each = [0] * num_subtypes
patient_types = {}
for line in ftype:
   words = line.strip().split()
   patient_types[words[0]] = words[1]
   subtype_index = subtypes.index(words[1])
   count_types_each[subtype_index] += 1

count_types_each_except_this = [0] * num_subtypes   
for i in range(num_subtypes):
   count_types_each_except_this = sum(count_types_each) - count_types_each[i]



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
   
for component in fhcomponents:
   genes = component.strip().split()
   patient_list_union = []
   for gene in genes: 
     patient_list_union.extend(dict_genes_to_patients[gene])
   patient_list_union = set(patient_list_union)
   count_subtypes_vec = count_subtype(patient_list_union)  
   for i in range(len(subtypes)):
      entry_11 = count_subtypes_vec[i]
      entry_12 = count_types_each[i] - count_subtypes_vec[i]
      entry_21 = sum(count_subtypes_vec) - count_subtypes_vec[i]
      entry_22 = count_types_each_except_this[i] - entry_21
      oddsratio, pvalue = stats.fisher_exact([[entry_11, entry_12], [entry_21, entry_22]])
      
      
      
      
      
      
      