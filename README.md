# MEXCOWalk: Mutual Exclusion and Coverage Based Random Walk to Identify Cancer Modules

### This is the original repository for MEXCOwalk paper


**Installing the dependencies**

```
pip install -r requirements.txt
```

## **Input**

### Graph Construction

There are two files for the graph construction: the network edges file and the nodes to index mappings
both files are located in data folder

1. The PPI network file:

the file is located at data/hint_edge_file.txt

```
node_i_Index node_j_Index weight
0 1 1
0 2 1
```

2. The nodes to index file mapping:

the file is located at data/hint_index_file.txt
```
index GeneName
1 A1BG
2 A1CF
```

### Mutation Data

The files needed for calculating the mutex and coverage scores are processed from src_snvs.tsv and src_cnas.tsv files.

Note: the files are tab separated

The files are:

1. Gene and its list of patients where that gene is mutated:
the file is located at data/genes_vs_patients_indices_gen_paran.txt

```
RNF14	109	241	...
UBE2Q1	581	837	...
....
```
the patients are referred to with their ids.


2. The patients name to patient id mapping:

the file is located at data/patients_to_indices_gen.txt

```
TCGA-06-0649	1
TCGA-BA-5149	2
TCGA-DK-A1A6	3
```
the ids starts from 1

3. The gene and its corresponding frequency:

This file contains the mutation frequncies, which are assigned as heats during the random walk.

the file is located at  data/pan12gene2freq.txt. the file was downloaded from [here](https://github.com/raphael-group/hotnet2/tree/master/paper/data/heats)

```
A1BG 0.00353697749196
A2M 0.0128617363344
A4GALT 0.00064308681672

```

### Cosmic genes

The cosmic genes are loaded from Census_allTue_May_23_12-08-15_2017.tsv file, which is located in data folder.



## **Run**

There are two bash scripts, one for running the modules discovery algorithm, the other for evaluation.

For more details on the execution parameters please refer to the bash files.

1. Modules discovery:

```
sh execute_all.sh
```
2. Evaluation

```
sh evaluate_all.sh
```

## **Knearest**


## **CoxPh analysis**
