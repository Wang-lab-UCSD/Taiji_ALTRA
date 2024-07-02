specificTFs_cell <- function(x,p,d,df,group,tag,ps=0,test="t.test",group_name="preClust"){
  
  n <- group[group[[group_name]]==x,"id"]
  n2 <- which(colnames(df) %in% n)
  print(length(n2))
    
  if (test=="t.test"){
      temp <- unpairedtTest(df, n2, p, d, ps=ps)
  } else if (test == "wilcoxon"){
      temp <- wilcoxonTest(df, n2, p, d, ps=ps)
  }
  Tcell <- rownames(temp)
  if (length(Tcell)==0){
    return()
  }
  ## write result
  fl <- paste0(x,"_p",p,"_d",d,"_gene",length(Tcell),"_",test,"_",tag,".txt")
  writeLines(Tcell, fl)
  
  fl2 <- paste0(x,"_p",p,"_d",d,"_gene",length(Tcell),"_",test,"_",tag,"_table.tsv")
  write.csv(temp, fl2)  
  return(Tcell)
}
