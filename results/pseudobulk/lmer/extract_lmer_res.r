files = grep(".Rdata", dir(), val=T)  
res = list() 
for( file in files){ 
    load(file)
    res[[file]] = colSums(fm.mat2[,grep("fdr", colnames(fm.mat2))] < 0.05, na.rm=T)  

}

temp = do.call(rbind, res) 
rownames(temp) = gsub("mod_ageDE_|sumcounts.|.Rdata", "", rownames(temp) ) 
colnames(temp) = c("sex+age", "age", "sex", "fem-age", "male-age")

write.table(temp, file="lmer_summary_results.txt")

