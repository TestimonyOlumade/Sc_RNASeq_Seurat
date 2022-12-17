## RUNNING DOUBLET FINDER ON ScRNA-SEQ DATA

# install DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

## Intro notes

# DoubletFinder needs 3 parameters to make doublet predictions
# - pN - nukber of artificial doubles. Default is 0.25 or 25%
# - pK - neighborhood size - to compute number of artificial nearest neighbors
# Exp - number of expected real doublets - can be gotten from user guide of reagent kits from 10x genomics

# you can use the no of cells loaded and no of cells recovered to estimate the 'Exp' value

## How the DoubletFinder algorithm works

# 1. Simulate doublets - averaging gene expression profiles of random pairs of cells,
# and introduces some proportion of artificial doublets into the dataset (pN)
# then artificial doublets are merged with real dataset using seurat scRNASeq pipeline 
# and standard pre-processing steps are followed

# 2. PCA - dimensionality reduction is performed on the merged real-artificial data using PCA
# this will produce a low-dimensional space describing similarity between real and artificial cells

# 3. Define the neighbors - the algorithm then detects the nearest neighbors for each cell in the PCA space
 
# 4. Calculate pANN threshold - this info to compute the proportion of artificial nearest neighbors (pANN) in each cell.
# # the performance of DoubletFinder is highly dependent on the pK parameter 

# 5. Remove doublets - real doublets are then predicted with top n pANN values
# where n is set to the total number of expected doublets


## Best practices

# not applied on aggregate or merged data
# run on distinct samples separately
# input should be cleared of low-quality cells
# clusters with low RNA UMIs, high % MT reads, and uninformative marker genes should be removed

------
## STEP 1 - read the raw feature and cell matrix data into R

# create counts matrix
cts <- ReadMtx(mtx = 'path_to_matrix_file.mtx.gz',
               features = 'path_to_features_file.tsv.gz',
               cells = 'path_to_barcodes_file.tsv.gz')

# to view the counts sparse matrix, looking at the first 10 rows and columns
cts[1:10,1:10]

# create Seurat object from the counts sparse matrix

pbmc.seurat <- CreateSeuratObject(counts = cts)

str(pbmc.seurat)

# QC and Filtering
# explore QC as with scRNA-Seq analysis using seurat

# estimate the % MT reads in cells

pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')

# this creates a column 'mitoPercent' in the seurat object

# then filter Seurat objects, keeping only desired cells for downstream analyses

pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

# to view the subset of Seurat object

pbmc.seurat.filtered

--------
## STEP 2: run standard pre-process standard workflow
  
# normalize data 

pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)

# find variable features

pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)

# scale data

pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)

# run linear dimensionality reduction

pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)

# find dimensionality of dataset 

ElbowPlot(pbmc.seurat.filtered)

# view dimensionality in ElbowPlot to see how variation is captured across the PCs
#then find neighbors using the first 20 dimensions, from ElbowPlot

pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)

# find clusters

pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)

# run UMAP

pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)


-------
## STEP 3: pK Identification (no ground-truth) --------------------

# finding the optimal pK value is paramount; it determines how well DoubletFinder works
# using the no ground-truth approach; when cell multiplexing data is not available 

# now introduce artificial doublets in various proportions and merge them with real dataset
# pre-process the dataset
# calculate the proportion of artificial nearest neighbors for varying neighborhood sizes

sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)

sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)

bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

# it finds the optimal pK values by computing the mean variance normalized by modality coefficient
# the result is a list of proportion of artificial nearest neighbors for varying combinations of pN and pK
# the pK value corresponding to the highest BCmetric value is the optimal pK value

# which can then be viewed using ggplot

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# now store the optimal pK value to a 'pK' variable

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 

# this gets returned as a list 
# so choose the first value on the list; convert to a numeric value to save the optimal pK value
optimal_pK <- as.numeric(as.character(pK[[1]]))


--------
## STEP 4: Homotypic Doublet Proportion Estimate ----------------------

# this can be estimated from user guide of reagent kits (no of cells loaded : no of cells recovered)
# this is done using user-provided annotations - which are cell clusters

annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters

# now model the proportion of homotypic doublets using the annotations

homotypic.prop <- modelHomotypic(annotations)  

# assuming a 7.6% doublet formation rate (from no of cells loaded and recovered)
# multiply that by the total number of cells in seurat object

nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  

# this returns the expected number of doublets
# now adjust the value for homotypic doublets

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# new expected doubled to be formed after adjustment for homotypic doublets


--------
## STEP 5: run doubletFinder 

pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = optimal_pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)

# if you have generated a pANN value from a previous run, 'pANN' in the code above can be set to TRUE
# if sct transform is done on the data during pre-processing, then 'sct' should be set to TRUE


--------
## STEP 6: visualize doublets

# the doublet predictions are stored in the metadata
  
view(pbmc.seurat.filtered@meta.data)

# a column is added as predictions for each cell - whether doublets or singlets
# the name of the column is then used to group the object (- use nrows)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")

# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)

# this returns the number of doublets to be filtered out


-------
## STEP 7: filter out and re-cluster the singlets
  
singlets <- subset(pbmc.seurat.filtered, subset = DF.classifications_0.25_0.21_691 == "Singlet")


singlets <- NormalizeData(object = singlets)

singlets <- FindVariableFeatures(object = singlets)

singlets <- ScaleData(object = singlets)

singlets<- RunPCA(object = singlets)

ElbowPlot(singlets)

singlets<- FindNeighbors(object = singlets, dims = 1:20)

singlets<- FindClusters(object = singlets)

singlets <- RunUMAP(object = singlets, dims = 1:20)





