bascetRoot <- getwd()

if(TRUE){
  
  setwd("~/github/zorn")
  
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R")
  source("R/zorn_aggr.R")
  source("R/count_kmer.R")
  source("R/countsketch.R")
  source("R/refgenome.R")
  source("R/kraken.R")
  source("R/container.R")
  source("R/ext_tools.R")
  
} else {
  library(Zorn)
}

################################################################################
################## Preprocessing with Bascet/Zorn ##############################
################################################################################

bascet_instance.default <- getBascetSingularityImage(store_at="~/mystore/")
bascet_runner.default <- SlurmRunner(account="hpc2n2025-074", ncpu="10")

bascetRoot <- file.path("/husky/henriksson/atrandi/250424_rnaseq6_miseq")  

setwd(bascetRoot)

###################################################
################################################### debarcode
###################################################

rawmeta_dir <- readLines("rawdata.txt")
rawmeta <- DetectRawFileMeta(rawmeta_dir)

barcode_error <- NULL
if(file.exists("barcode_error.txt")) {
  barcode_error <- readLines("barcode_error.txt")
}

system("echo START BascetGetRaw >> time.txt; echo `date +%s` >> time.txt")
BascetGetRaw(
  bascetRoot,
  rawmeta,
  chemistry="atrandi_rnaseq"
)
system("echo END BascetGetRaw >> time.txt; echo `date +%s` >> time.txt")

###################################################
################################################### shardify
###################################################

### Decide cells to include
h <- ReadHistogram(bascetRoot,"debarcoded")
includeCells <- h$cellid[h$count>10]       ########### 10 for miseq
#includeCells <- h$cellid[h$count>100]       ########### 10 for miseq
length(includeCells)

### Shardify i.e. divide into multiple sets of files for parallel processing.
# this command will spawn many processes, as it does random I/O on the input files
system("echo START BascetShardify >> time.txt; echo `date +%s` >> time.txt")
BascetShardify(
  bascetRoot,
  includeCells = includeCells,
  num_output_shards = 1, ########### up to 20
  runner=SlurmRunner(bascet_runner.default, ncpu="4")  #not much CPU needed. increased for memory demands
)
system("echo END BascetShardify >> time.txt; echo `date +%s` >> time.txt")



### Get reads in fastq format for BWA
system("echo START asfq >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  inputName="filtered",
  outputName="asfq",
  out_format="R1.fq.gz",
  runner=SlurmRunner(bascet_runner.default, ncpu="4")  #in and out processes
)
system("echo END asfq >> time.txt; echo `date +%s` >> time.txt")



### Get reads in fastq format for BWA
system("echo START fastp >> time.txt; echo `date +%s` >> time.txt")
BascetRunFASTP(
  bascetRoot,
  numLocalThreads=10,
  inputName="asfq",
  outputName="fastp"
)
system("echo END fastp >> time.txt; echo `date +%s` >> time.txt")



system("echo START totirp >> time.txt; echo `date +%s` >> time.txt")
BascetMapTransform(
  bascetRoot,
  "fastp",
  "new_filtered",
  out_format="tirp.gz"
)
system("echo END totirp >> time.txt; echo `date +%s` >> time.txt")




################################################################################
################## Alignment ###################################################
################################################################################


### Perform alignment
system("echo START BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")
BascetAlignToReference(
  bascetRoot,
  useReference="/home/m/mahogny/mystore/atrandi/bwa_ref/yersinia/all.fa",
  numLocalThreads=10
)
system("echo END BascetAlignToReference >> time.txt; echo `date +%s` >> time.txt")



### Generate fragments BED file suited for quantifying reads/chromosome using Signac later -- this is a wrapper for mapshard
system("echo START BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")
BascetBam2Fragments(
  bascetRoot,
  runner=SlurmRunner(bascet_runner, ncpu="2")
)
system("echo END BascetBam2Fragments >> time.txt; echo `date +%s` >> time.txt")




