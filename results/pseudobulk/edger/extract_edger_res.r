files = grep(".Rdata", dir(), val=T)  
res = list() 
for( file in files){ 
    load(file)

    temp1 = top.table.list[[1]] 
    temp2 = top.table.list[[2]] 
    temp3 = top.table.list[[3]] 
    temp4 = top.table.list[[4]] 
    temp5 = top.table.list[[5]] 
    temp6 = top.table.list[[6]] 
    temp7 = top.table.list[[7]] 
    
    o = order(rownames(temp1))
    temp1 = temp1[o,]

    o = order(rownames(temp2))
    temp2 = temp2[o,]

    o = order(rownames(temp3))
    temp3 = temp3[o,]

    o = order(rownames(temp4))
    temp4 = temp4[o,]

    o = order(rownames(temp5))
    temp5 = temp5[o,]

    o = order(rownames(temp6))
    temp6 = temp6[o,]
    o = order(rownames(temp7))
    temp7 = temp7[o,]

    temp_logFC = abs(temp1$logFC) > 1 
 
    res[[file]] = c(sum(temp1$adj.P.Val < 0.05 & temp_logFC  , na.rm=T),
    sum(temp2$adj.P.Val < 0.05 & temp_logFC , na.rm=T),
    sum(temp3$adj.P.Val < 0.05 & temp_logFC, na.rm=T),
    sum(temp4$adj.P.Val < 0.05 & temp_logFC, na.rm=T),
    sum(temp5$adj.P.Val < 0.05 & temp_logFC, na.rm=T),
    sum(temp6$adj.P.Val < 0.05 & temp_logFC, na.rm=T),
    sum(temp7$adj.P.Val < 0.05 & temp_logFC, na.rm=T))   

}

 
temp = do.call(rbind, res) 
rownames(temp) = gsub("limma_|sumcounts.|.Rdata", "", rownames(temp) ) 
colnames(temp) = c("sex", "age", "pool", "sex+age+pool_trend", "sex+age+pool", "sex+age-pool", "sex+age+-pool_trend" )
 
write.table(temp, file="edger_summary_results.txt")

