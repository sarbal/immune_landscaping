## Extract data from larger Seurat object. Example code below. 
## For "raw" analyses - used OBJ@assays$RNA@counts, for normalised analyses, used OBJ@assays$SCT@counts
## Note, generated these prior to running MAST and stored as the large Seurat object took up lots of memory and failed to run on cluster
##


##    f1 = obj@meta.data[[opt$cell_level]] == celltypei & !is.na(obj@meta.data[[opt$cell_level]])
##    pbmc = obj[,f1]      
##    counts =  pbmc@assays$RNA@counts
##    bin_counts =  counts > 0
##    rm(counts)
##    N =  dim(bin_counts)[2]
##    n_cells = rowSums(bin_counts) 
##    counts_freq = n_cells/N 
##    DefaultAssay(pbmc) <- "RNA"
##    sca <- as(as.SingleCellExperiment(pbmc), 'SingleCellAssay')
##    ## MAST works on log2-normalized count matrix (see equivalent matrices of counts/logcounts in Seurat object vs. SingelCellAssay object)
##    assayNames(sca) #which assays in sca --> counts and logcounts
##    sca <- sca[counts_freq>0,] #returns the frequency of expression, i.e., the proportion of non-zero values in sc
##    colData(sca)$wellKey <- colnames(sca)
##    rowData(sca)$primerid <- rownames(rowData(sca))
##    colData(sca) <- droplevels(colData(sca))
##    colData(sca)$age_scaled <- scale(colData(sca)$age) # rescale continous variable (Age)
##    cdr2 <- colSums(assay(sca)>0) #assay(sca) is working on assays(sca)$logcounts
##    colData(sca)$cngeneson <- scale(cdr2)
##    save(sca , file=paste0(celltypei, "_scaraw.Rdata") ) 


## MAST code
## Note, sample code below was run in script where celltypei was specified and other variables pre-loaded 

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
 

 