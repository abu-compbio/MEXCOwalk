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
```
node_i_Index node_j_Index weight
0 1 1
0 2 1
```

2. The nodes to index file mapping:
```
index GeneName
1 A1BG
2 A1CF
```

### Mutations Data

The files needed for calculating the mutex and coverage scores are processed from src_snvs.tsv and src_cnas.tsv files

The files are:

NOTE: the files are tab separated

1. The list of patients whose gene ~~i~~ is mutated:
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

### Cosmic genes

The comic genes are loaded from Census_allTue_May_23_12-08-15_2017.tsv file located in data



## **Run**

There are two bash scripts, one for running the modules discovery algorithm, the other for evaluation
For more details on the execution parameters please refer to the bash files.

1. Modules discovery:

```
sh execute_all.sh
```
2. Evaluation

```
sh evaluate_all.sh
```
