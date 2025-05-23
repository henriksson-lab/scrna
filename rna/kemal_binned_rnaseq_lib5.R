
#could not find function "BascetGetRaw"

library(Signac)
library(Seurat)
library(ggplot2)

### TODO recompute rnaseq3 too

#for(curlib in c(2:3)){
#  curlib 
  
  ################################################################################
  ################## Preprocessing with Bascet/Zorn ##############################
  ################################################################################
  
  inst_bascet <- LocalRunner(direct = TRUE, show_script=TRUE)
  bascetRoot <- file.path("/husky/henriksson/atrandi/rnaseq5/frag")  ###change
  fastqdir <- file.path("/husky/fromsequencer/250311_failed_saliva/joram_rnaseq/raw")  
  rawmeta <- DetectRawFileMeta(fastqdir, verbose = TRUE)
  # 
  # print(bascetRoot)
  # print(rawmeta)
  # 
  # bascet_instance.default  #temp dir is here
  # 
  # 
  # ### Debarcode the reads, then sort them.
  # BascetGetRaw(
  #   bascetRoot,
  #   rawmeta,
  #   runner=inst
  # )
  # 
  # ### Decide cells to include
  # h <- ReadHistogram(bascetRoot,"debarcoded")
  # PlotHistogram(h)
  # includeCells <- h$cellid[h$count>100]   ### for rnaseq miseq #3.1
  # length(includeCells)
  # 
  # ### Shardify i.e. divide into multiple sets of files for parallel processing
  # BascetShardify(
  #   bascetRoot,
  #   #includeCells = includeCells, ## include all cells! filter later   TODO it is really expensive to search for individual cells. need other mode for merging
  #   runner = inst
  # )
  # 
  # 
  # 
  # ### Note: shardification can partially be skipped for RNAseq if you know what you are doing.
  # ### but this is poorly tested, so doing it to be sure!
  # 
  # 
  # 
  # 
  # 
  # 
  # ################################################################################
  # ################## Reference-based mapping to get "ground truth" ###############
  # ################################################################################
  # 
  # 
  # ### Get reads in fastq format for BWA  --- will not be needed later as the conversion will be made on the fly
  # BascetMapTransform(
  #   bascetRoot, 
  #   "filtered", 
  #   "asfq",
  #   out_format="fq.gz",   ## TODO we need two fq as out!! ideally at least. or if R1.fq.gz => write two of them. otherwise gather?
  #   runner=inst
  # )
  # 
  # ### Perform alignment -- internally wraps mapshard
  # BascetAlignToReference(  #### does not run in parallel?? TODO bug. use runjob
  #   bascetRoot,
  #   useReference="/husky/fromsequencer/241210_joram_rnaseq/ref/all.fa",
  #   numLocalThreads=10
  # )
  # 
  # #"/husky/fromsequencer/241210_joram_rnaseq/ref/all.fa"
  # 
  # ### Generate fragments BED file suited for quantifying reads/chromosome using Signac later
  # BascetBam2Fragments(
  #   bascetRoot,
  #   runner=inst
  # )
  
  
  ################################################################################
  ################ Use Signac/Seurat as a feature counter ########################
  ################################################################################
  
  
  
  
  ####### Find coordinates of each gene for our organism of choice
  gff <- rtracklayer::readGFF("/husky/fromsequencer/241210_joram_rnaseq/ref/all.gtf")
  #gff <- rtracklayer::readGFF("/husky/fromsequencer/241210_joram_rnaseq/ref/all.gff3")
  gff <- gff[gff$type == "transcript",]
#  gff$Name <- gff$gene_id
#  gff$Name[!is.na(grange_gene$gene_name)] <- gff$gene_name[!is.na(grange_gene$gene_name)]
  
  ### all ID??
  grange_gene <- GenomicRanges::makeGRangesFromDataFrame(gff)
  #grange_gene$type <- gff$type
  #grange_gene$gene_biotype <- gff$gene_biotype
  grange_gene$Name <- gff$gene_id
  grange_gene$Name[!is.na(gff$gene_name)] <- gff$gene_name[!is.na(gff$gene_name)]
  
  curlib <- "1"
  
  ####### Perform counting
  adata_f <- FragmentsToSignac(file.path(bascetRoot,paste0("fragments.",curlib,".tsv.gz")))  ## TODO, support multiple fragment files
  #adata[["RNA"]] <- CountGrangeFeatures(adata, grange_gene) ## this fails for no good reason; replacement has 8799 rows, data has 2
  adata <- CreateSeuratObject(CountGrangeFeatures(adata_f, grange_gene))
  
  ####### 
  ####### Divide libraries
  ####### 
  
  adata$well_r1 <- stringr::str_split_i(colnames(adata),"_",2)
  
  adata$libname <- curlib
  
  saveRDS(adata, file.path(bascetRoot,"adata2.RDS"))
  
  
  
  
  ################################################################################
  ################ Use Signac/Seurat as a bin counter ########################
  ################################################################################
  
  
  
  
  ####### Find coordinates of each gene for our organism of choice
  gff <- rtracklayer::readGFF("/husky/fromsequencer/241210_joram_rnaseq/ref/all.gtf")
  chrom_len <- sqldf::sqldf("select seqid, max(end) as len from gff group by seqid")
  bin_edge <- round(seq(from=1, to=4690288, length.out=40))
  
  bins1 <- data.frame(seqid="NZ_CP009792.1", start=bin_edge[1:(length(bin_edge)-1)], end=bin_edge[2:length(bin_edge)])
  bins2 <- data.frame(seqid="NC_006153.2", start=1, end=68141)
  
  gff <- rbind(bins1,bins2)
  
  grange_gene <- GenomicRanges::makeGRangesFromDataFrame(gff)
  #grange_gene$type <- gff$type
  #grange_gene$gene_biotype <- gff$gene_biotype
  grange_gene$Name <- paste0(gff$seqid,"_",1:nrow(gff))
  #grange_gene$Name[!is.na(gff$gene_name)] <- gff$gene_name[!is.na(gff$gene_name)]
  
  curlib <- "1"
  
  ####### Perform counting
  adata_f <- FragmentsToSignac(file.path(bascetRoot,paste0("fragments.",curlib,".tsv.gz")))  ## TODO, support multiple fragment files
  #adata[["RNA"]] <- CountGrangeFeatures(adata, grange_gene) ## this fails for no good reason; replacement has 8799 rows, data has 2
  adata <- CreateSeuratObject(CountGrangeFeatures(adata_f, grange_gene))
  
  ####### 
  ####### Divide libraries
  ####### 
  
  adata$well_r1 <- stringr::str_split_i(colnames(adata),"_",2)
  
  adata$libname <- curlib
  
  saveRDS(adata, file.path(bascetRoot,"binned_adata2.RDS"))
  

####### TODO generate 100kb bins instead
  

################################################################################ 
########################### Analyze merged count files ######################### 
################################################################################ 



#bascetRoot <- file.path("/husky/henriksson/atrandi/rnaseq3", curlib)

### 3 aliquots; we can just concatenate these adata files

adata1 <- readRDS("/husky/henriksson/atrandi/rnaseq5/adata.RDS")
adata2 <- readRDS("/husky/henriksson/atrandi/rnaseq5/adata2.RDS")
adata3 <- readRDS("/husky/henriksson/atrandi/rnaseq5/adata.RDS")

#adata1$libname <- "1"
#adata2$libname <- "2"
#adata3$libname <- "3"

pbmc.combined <- merge(
  adata1, y = c(adata2, adata3)
, add.cell.ids = c("sub1", "sub2","sub3"), project = "all")

pbmc.combined <- adata2

##### Assign metadata to each cell from sample layout file
samplelayout <- read.csv("/husky/fromsequencer/250204_joram_rnaseq4_miseq/raw/sample_layout.csv", sep="\t")
rownames(samplelayout) <- samplelayout$well
adata$libname <- samplelayout[adata$well_r1,]$lib


#saveRDS(adata, file.path("/husky/henriksson/atrandi/rnaseq5/combined.RDS"))
#saveRDS(adata, file.path("/husky/henriksson/atrandi/rnaseq5/combined.RDS"))


################################################################################ 
########################### Kneeplot analysis ################################## 
################################################################################ 

### use cached data
#adata <- readRDS("/husky/henriksson/atrandi/rnaseq5/combined.RDS")
adata

####### 
####### Perform knee plot, all libs 
####### 

### All cells together  --- counts
df <- data.frame(cnt=sort(adata$nCount_RNA, decreasing = TRUE))
df$rank <- 1:nrow(df)
ggplot(df,aes(rank, cnt)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()

### All cells together --- number of genes
df <- data.frame(cnt=sort(adata$nFeature_RNA, decreasing = TRUE))
df$rank <- 1:nrow(df)
ggplot(df,aes(rank, cnt)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()

####### 
####### Perform knee plot, for each lib
####### 

all_df <- NULL
for(curlib in unique(adata$libname)){
  df <- data.frame(cnt=sort(adata$nCount_RNA[adata$libname==curlib], decreasing = TRUE), lib=curlib)
  df$rank <- 1:nrow(df)
  all_df <- rbind(all_df, df)
}
ggplot(all_df,aes(rank, cnt, color=lib)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()



####### 
####### Perform RNA-seq type analysis; see https://satijalab.org/seurat/articles/adata3k_tutorial.html
####### 

DefaultAssay(adata) <- "RNA"

#hist(log10(adata$nCount_RNA))
adata <- adata[,adata$nCount_RNA>100] #meaningless below this number of reads for certain!!
dim(adata)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adata), 20)
top10

# plot variable features with and without labels
if(FALSE){
  plot1 <- VariableFeaturePlot(adata)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
}


all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)

adata <- RunPCA(adata, features = rownames(adata)[rowSums(adata@assays$RNA@counts)<50000])#VariableFeatures(object = adata))
#adata <- RunPCA(adata, features = all.genes)#VariableFeatures(object = adata))
adata <- FindNeighbors(adata, dims = 1:20)
adata <- RunUMAP(adata, dims = 1:20)

adata <- FindClusters(adata, resolution = 0.2)

DimPlot(adata, reduction = "umap")

#DimPlot(adata, reduction = "umap", group.by = "libname")


FeaturePlot(adata, features=top10)

FeaturePlot(adata, features=rownames(adata)[1:10])
FeaturePlot(adata, features="NC-006153.2-40" )

######## Names of top genes
sort(rowSums(adata@assays$RNA@counts), decreasing = TRUE)[1:100]

df <- data.frame(
  cnt = sort(rowSums(adata@assays$RNA@counts), decreasing = TRUE)
)
df$index <- 1:nrow(df)
ggplot(df, aes(index, cnt)) + geom_point() + scale_y_log10()
#plot(sort(rowSums(adata@assays$RNA@counts), decreasing = TRUE))

df <- data.frame(
  cnt = rowSums(adata@assays$RNA@counts),
  gene = rownames(adata@assays$RNA@counts)
)
df <- df[order(df$cnt, decreasing = TRUE),]
df$rank <- 1:nrow(df)
ggplot(df, aes(rank, cnt)) + geom_point() + scale_y_log10()

df

df$gene[1:38]

FeaturePlot(adata, features=unique(names(sort(rowSums(adata@assays$RNA@counts), decreasing = TRUE)))[25:35])

########## remove top genes
bdata <- adata[!(rownames(adata) %in% df$gene[1:10]),]
bdata <- RunPCA(bdata)#, features = VariableFeatures(object = bdata))
bdata <- FindNeighbors(bdata, dims = 1:10)
bdata <- RunUMAP(bdata, dims = 1:10)

bdata <- FindClusters(bdata, resolution = 0.1)

DimPlot(bdata, reduction = "umap")
bdata


########## remove BZ22 genes
bdata <- adata[!str_detect(rownames(adata), "BZ22"),]
bdata
bdata <- RunPCA(bdata, features = rownames(bdata))
bdata <- FindNeighbors(bdata, dims = 1:50)
bdata <- RunUMAP(bdata, dims = 1:50)
bdata <- FindClusters(bdata, resolution = 0.1)
DimPlot(bdata, reduction = "umap")


######### compare with featurecounts results
df <- data.frame(
  cnt = sort(as.integer(read.table("/husky/henriksson/atrandi/rnaseq5/testal.cnt",sep="\t")$V7[-1]), decreasing = TRUE)
)
df <- df[order(df$cnt, decreasing = TRUE),,drop=FALSE]
df$rank <- 1:nrow(df)
ggplot(df, aes(rank, cnt)) + geom_point() + scale_y_log10()
  



if(FALSE){
  library(dplyr)
  pbmc.markers <- FindAllMarkers(bdata, only.pos = FALSE)

  df <- as.data.frame(pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.2) %>%
    dplyr::filter(pct.1 > 0.01))
  df$Name <- df$gene
  #df <- merge(df, gff[,c("Name","type")])
  df <- df[order(df$p_val),]
  df

}


FeaturePlot(adata, features=c("gene-YPTB-RS21540"))
FeaturePlot(adata, features=c("gene-BZ22-RS01385"))
FeaturePlot(adata, features=c("icd"))
FeaturePlot(adata, features=c("odhB"))



FeaturePlot(adata, features=df$gene[1:10])



FeaturePlot(adata, features=c("ssrA")) #transfer-messenger RNA
FeaturePlot(adata, features=c("BZ22-RS01385")) #another major driver. rRNA

FeaturePlot(adata, features=c("yqhD"))
FeaturePlot(adata, features=c("rnpB"))
FeaturePlot(adata, features=c("ompA"))
FeaturePlot(adata, features=c("gndA"))


FeaturePlot(adata, features=c("katG"))
FeaturePlot(adata, features=c("ssrA")) #drives a large part of the clustering. transfer-messenger RNA"

FeaturePlot(adata, features=c("BZ22-RS04865"))



####### 
####### QC for each library (1-3)
####### 

DimPlot(adata, reduction = "umap", group.by = "libname")  #lib2 enriched in ssrA cluster

DimPlot(adata, reduction = "umap", group.by = "well_r1")#, label = TRUE)



####### 
####### For each library, count features
####### 

#curlib <- "lib1"
alldf <- NULL
for(curlib in unique(adata$libname)){
  df <- as.data.frame(merge(
    data.frame(
      Name=rownames(adata),
      cnt=rowSums(adata[,adata$libname==curlib]@assays$RNA@counts),
      lib=curlib
    ),
    as.data.frame(grange_gene)[,c("Name","gene_biotype")]
  ))
  #as.data.frame(df)
  alldf <- rbind(
    alldf,
    sqldf::sqldf("select gene_biotype, sum(cnt) as cnt, lib from df group by gene_biotype")
  )
  #protein_coding is the relevant one  
}
sum_cnt_type <- reshape2::acast(data = alldf, gene_biotype~lib, value.var = "cnt")

sum_cnt_type_norm <- sum_cnt_type
for(i in 1:ncol(sum_cnt_type_norm)){
  sum_cnt_type_norm[,i] <- round(digits = 3,sum_cnt_type_norm[,i]/sum(sum_cnt_type_norm[,i]))
}

sum_cnt_type
sum_cnt_type_norm




####### 
####### Pileups
####### 

adata_f <- CreateSeuratObject(adata_f)
adata_f

adata_f$nCount_RNA <- adata$nCount_RNA
adata_f$libname <- adata$libname

CoveragePlot(
#  assay.scale = "separate",
  object = adata_f,
  region = "NZ_CP009792.1-60000-90000", 
  group.by = "libname"
)















###################
######## remove rRNA
###################




all.genes <- rownames(bdata)
bdata <- ScaleData(bdata, features = all.genes)
bdata <- RunPCA(bdata, features = VariableFeatures(object = bdata))
bdata <- FindNeighbors(bdata, dims = 1:10)
bdata <- RunUMAP(bdata, dims = 1:10)

bdata <- FindClusters(bdata, resolution = 0.2)

DimPlot(bdata, reduction = "umap")

DimPlot(bdata, reduction = "umap", group.by = "libname")

FeaturePlot(bdata, features=unique(names(sort(rowSums(bdata@assays$RNA@counts), decreasing = TRUE)))[1:5])



