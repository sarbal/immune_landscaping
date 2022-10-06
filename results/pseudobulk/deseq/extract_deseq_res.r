files = grep(".Rdata", dir(), val=T)  
res = list() 
  thrsh =  0.1 
  fdr = 0.05 
for( file in files){ 
    load(file)

   res[[file]]  =  sum(degs$padj2  < fdr & abs(degs$log2_fc) > thrsh ),
    sum(degs2$padj2 < fdr & abs(degs2$log2_fc) > thrsh ),
    sum(degs3$padj2 < fdr & abs(degs3$log2_fc) > thrsh ))


}

 
temp = do.call(rbind, res) 
rownames(temp) = gsub("DESeq_age_|sumcounts.|.Rdata", "", rownames(temp) ) 
colnames(temp) = c("sex+age_scaled", "reduced-sex", "contrast-sex")
 
write.table(temp, file="deseq_summary_results.txt")

