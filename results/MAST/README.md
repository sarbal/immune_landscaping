# Summary 

|                                                         | Number of cells (total) | Number of cells (male) | Number of cells (female) | Number of DE genes (sex) | Number of DE genes (age) |
|---------------------------------------------------------|-------------------------|------------------------|--------------------------|--------------------------|--------------------------|
|                                                         |                         |                        |                          |                          |                          |
| Cell-types   (l1)                                       |                         |                        |                          |                          |                          |
| B                                                       | 131968                  | 48632                  | 83336                    | 28                       | 2                        |
| CD4 T                                                   | 624800                  | 249404                 | 375396                   | 267                      | 4                        |
| CD8 T                                                   | 243186                  | 104407                 | 138779                   | 196                      | 57                       |
| DC                                                      | 6511                    | 3201                   | 3310                     | 12                       | 0                        |
| Mono                                                    | 54034                   | 25648                  | 28386                    | 36                       | 0                        |
| NK                                                      | 183022                  | 84157                  | 98865                    | 134                      | 3                        |
| other T                                                 | 20355                   | 8105                   | 12250                    | 19                       | 0                        |
| other                                                   | 3892                    | 1713                   | 2179                     | 6                        | 0                        |
|                                                         |                         |                        |                          |                          |                          |
| Cell-types   (l2)                                       |                         |                        |                          |                          |                          |
| B intermediate                                          | 33194                   | 12440                  | 20754                    | 15                       | 0                        |
| B memory                                                | 23708                   | 8582                   | 15126                    | 12                       | 0                        |
| B naive                                                 | 71374                   | 25965                  | 45409                    | 18                       | 0                        |
| CD14 Mono                                               | 38908                   | 18732                  | 20176                    | 27                       | 0                        |
| CD16 Mono                                               | 15341                   | 7011                   | 8330                     | 19                       | 0                        |
| CD4 CTL                                                 | 11680                   | 5049                   | 6631                     | 11                       | 1                        |
| CD4 Naive                                               | 303555                  | 114708                 | 188847                   | 157                      | 8                        |
| CD4 Proliferating                                       | 556                     | 251                    | 305                      | 7                        | 2                        |
| CD4 TCM                                                 | 251068                  | 106156                 | 144912                   | 119                      | 3                        |
| CD4 TEM                                                 | 24266                   | 9665                   | 14601                    | 17                       | 1                        |
| CD8 Naive                                               | 57914                   | 22248                  | 35666                    | 48                       | 0                        |
| CD8 TCM                                                 | 10724                   | 4758                   | 5966                     | 15                       | 1                        |
| CD8 TEM                                                 | 178401                  | 79043                  | 99358                    | 64                       | 2                        |
| DC**                                                    | 6596                    | 3227                   | 3369                     | 12                       | 0                        |
| Eryth                                                   | 580                     | 232                    | 348                      | 5                        | 2                        |
| HSPC                                                    | 1766                    | 864                    | 902                      | 5                        | 0                        |
| MAIT                                                    | 12286                   | 4511                   | 7775                     | 10                       | 0                        |
| NK Proliferating                                        | 3215                    | 1579                   | 1636                     | 10                       | 1                        |
| NK_CD56bright                                           | 8259                    | 3112                   | 5147                     | 13                       | 0                        |
| NK                                                      | 171757                  | 79541                  | 92216                    | 117                      | 4                        |
| Plasmablast                                             | 3602                    | 1618                   | 1984                     | 13                       | 0                        |
| Treg                                                    | 28843                   | 11502                  | 17341                    | 23                       | 2                        |
|                                                         |                         |                        |                          |                          |                          |


 
| ** DC cells were combined from these classifications: |      |
|-------------------------------------------------------|------|
| cDC1                                                  | 112  |
| cDC2                                                  | 4330 |
| ASDC                                                  | 210  |
| pDC                                                   | 1944 |
|                                                       |      |
| Cell   types not assessed:                            |      |
| ILC                                                   | 315  |
| Platelet                                              | 1773 |
| gdT                                                   | 4883 |
| dnT                                                   | 2939 |
| CD8 Proliferating                                     | 256  |


# MAST code/pseudocode
First extracted data from larger Seurat object. Example code below. 
- For "raw" analyses - used OBJ@assays$RNA@counts, for normalised analyses, used OBJ@assays$SCT@counts
- Note, generated these prior to running MAST and stored as the large Seurat object took up lots of memory and failed to run on cluster
- This data is not shared (large, private)

```
f1 = obj@meta.data[[opt$cell_level]] == celltypei & !is.na(obj@meta.data[[opt$cell_level]])
pbmc = obj[,f1]      
counts =  pbmc@assays$RNA@counts
bin_counts =  counts > 0
rm(counts)
N =  dim(bin_counts)[2]
n_cells = rowSums(bin_counts) 
counts_freq = n_cells/N 
DefaultAssay(pbmc) <- "RNA"
sca <- as(as.SingleCellExperiment(pbmc), 'SingleCellAssay')
## MAST works on log2-normalized count matrix (see equivalent matrices of counts/logcounts in Seurat object vs. SingelCellAssay object)
assayNames(sca) #which assays in sca --> counts and logcounts
sca <- sca[counts_freq>0,] #returns the frequency of expression, i.e., the proportion of non-zero values in sc
colData(sca)$wellKey <- colnames(sca)
rowData(sca)$primerid <- rownames(rowData(sca))
colData(sca) <- droplevels(colData(sca))
colData(sca)$age_scaled <- scale(colData(sca)$age) # rescale continous variable (Age)
cdr2 <- colSums(assay(sca)>0) #assay(sca) is working on assays(sca)$logcounts
colData(sca)$cngeneson <- scale(cdr2)
save(sca , file=paste0(celltypei, "_scaraw.Rdata") ) 
```

## MAST analysis
Note, sample code below was run in script where celltypei was specified and other variables pre-loaded.
- Shared the "fcHurdle" files in folders. 

## Step 1 - load data 
    load(file=paste0(celltypei, "_scaraw.Rdata") ) 

## Step 2 - setup remaining variables 
    opt[['cell_type']] = celltypei
    opt[['cell_level']] = "predicted.celltype.l1" ## For Azimuth L1 
    opt[['cell_level']] = "predicted.celltype.l2" ## For Azimuth L2 

    cell_type1 <- opt$cell_type
    cell_level <- opt$cell_level

## Step 3 - setup sca object for analysis, this includes removing empty cells/samples and low frequency genes 
    counts =  assays(sca)$counts
    bin_counts =  counts > 0
    rm(counts)
    N =  dim(bin_counts)[2]
    n_cells = rowSums(bin_counts) 
    counts_freq = n_cells/N
    rm(bin_counts)

    # Apply function to 'sca' object testing contrast.vars (from opt$covs, in this case 'Age' and 'Age_scaled')
    contrast.vars <- unique(covs.df[covs.df$type=='contrast',]$covariate)
    fixed_effects <- unique(covs.df[covs.df$type=='fixed',]$covariate)
    random_effects <- unique(covs.df[covs.df$type%in%c('random','batch'),]$covariate)


    # counts_freq = counts_freq[counts_freq>0]
    contrast = contrast.vars[1]
    freq_expressed <- 0.1
    print(paste0('Setting a minimum gene expression threshold based on genes that are found in at least ', freq_expressed, ' of the cells (proportion of non-zero cells)'))
    expressed_genes <- counts_freq > freq_expressed  
    sca_object <- sca[expressed_genes,] 


    print(paste0('Testing: ', contrast))
    random_effects.fmla <- paste(paste0('(1|',random_effects,')'),collapse='+') #try other nomenclature (=option 1, default)
    contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
    zlm_vars <- paste0('~',paste(c(contrast_fixed.fmla,random_effects.fmla), collapse='+'))
    zlm_formula <- as.formula(zlm_vars)
     
## Step 4 - save progress 
    zlmfile=paste0("zlm_freq01_", opt[['cell_type']] ,flag, ".Rdata")


## Step 5 - run conditional analysis. Due to failing scripts this is run once (and saved) or loaded for downstream contrasts 
    if( !file.exists ( zlmfile   ) )  {  
        print(paste0('Trying to mitigate model failing to converge for some genes by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))'))
        zlmCond <- zlm(zlm_formula, sca_object, method='glmer', ebayes=FALSE, fitArgsD=list(nAGQ=0))
        save(zlmCond, file=zlmfile) 
    } else { 
        load(zlmfile)
    }

## Step 6 - age contrast 
    summaryfile=paste0("summaryDT_batch_age_",opt[['cell_type']] ,flag,".rds")
    if( !file.exists ( summaryfile   ) )  {  
        summaryCond <- summary(zlmCond, doLRT=contrast, fitArgsD=list(nAGQ=0))
        summaryDt <- summaryCond$datatable
        saveRDS(summaryDt, file=summaryfile) 

        contrast <- 'sex'
        fcHurdle = summaryDt_sex[ summaryDt_sex$contrast == contrast & summaryDt_sex$component=='H', c(1,4)]
        fcHurdle2 = summaryDt_sex[ summaryDt_sex$contrast == contrast & summaryDt_sex$component=='logFC', c(7,5,6)]
        fcHurdle_sex = cbind(fcHurdle, fcHurdle2, p.adjust(fcHurdle[,2]) )
        colnames(fcHurdle_sex)[6] = "fdr"
        save(fcHurdle_sex, file="fcHurdle_sex",opt[['cell_type']] ,flag,".Rdata") 
    } 

## Step 7 - sex contrast 
    summaryfile=paste0("summaryDT_batch_sex_",opt[['cell_type']] ,flag,".rds")
    if( !file.exists ( summaryfile   ) )  {  
        summaryCond <- summary(zlmCond, doLRT="sex", fitArgsD=list(nAGQ=0))
        summaryDt <- summaryCond$datatable
        saveRDS(summaryDt, file=summaryfile) 

        contrast <- 'age'
        fcHurdle = summaryDt_f[ summaryDt_f$contrast == contrast & summaryDt_f$component=='H', c(1,4)]
        fcHurdle2 = summaryDt_f[ summaryDt_f$contrast == contrast & summaryDt_f$component=='logFC', c(7,5,6)]
        fcHurdle_age = cbind(fcHurdle, fcHurdle2, p.adjust(fcHurdle[,2]) )
        colnames(fcHurdle_age)[6] = "fdr"
        save(fcHurdle_age, file="fcHurdle_age",opt[['cell_type']] ,flag,".Rdata") 

    } 


 
