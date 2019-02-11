knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)
devtools::load_all(".")
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(TCGAbiolinks)
library(DESeq2)

dir.create("../data/TCGA/"), showWarnings = FALSE)

cancers = c('BRCA','BLCA','CRC', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC', 'COAD')


for (cancer in cancers){

  CancerProject <- sprintf("TCGA-%s", cancer)
  DataDirectory <- paste0("../data/TCGA/",gsub("-","_",CancerProject))
  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")


  query <- GDCquery(project = CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts")

  samplesDown <- getResults(query,cols=c("cases"))
  if (cancer == 'LAML') {
    smTP = 'TB'
  } else{
    smTP = 'TP'
  }
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                    typesample = smTP)



  queryDown <- GDCquery(project = CancerProject,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts",
                        barcode = c(dataSmTP))

  GDCdownload(query = queryDown,
              directory = DataDirectory)


  dataPrep <- GDCprepare(query = queryDown,
                         save = TRUE,
                         directory =  DataDirectory,
                         save.filename = FileNameData)

  #load(FileNameData)
  data_<- assay(dataPrep)

  ids <- data_[1,]
  condition <- rep('treated', length(ids))
  type <- rep('SOMETYPE', length(ids))

  coldata <- cbind(ids, condition, type)
  ss <- coldata[, c(2,3)]
  dds <- DESeqDataSetFromMatrix(countData = data_,
                                colData = ss,
                                design = ~ 1)

  dds <- estimateSizeFactors(dds)
  rld <- vst(dds, blind=TRUE)
  res <- assay(rld)

  SaveDirectory <- paste0("../data/TCGA/",cancer)
  FileNameSave <- paste0(SaveDirectory, "_","Normalized",".csv")
  write.csv(as.data.frame(res), file=FileNameSave)

}
