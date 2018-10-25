files = grep("*virome_counts_unique_noTanRep_final.txt",list.files(),value=TRUE)
cnts = lapply(files,function(i){
  read.csv(i,sep="\t",header=FALSE,stringsAsFactors=FALSE,row.names=1)
})
viruses = unique(unlist(lapply(cnts,function(i)rownames(i))))
samples = paste(sapply(strsplit(files,"\\_virome"),function(i)i[1]),sep="_")
counts = array(NA,dim=c(length(viruses),length(samples)))
colnames(counts) = samples
rownames(counts) = viruses
for(i in seq(cnts)){
  counts[rownames(cnts[[i]]),i] = cnts[[i]][,2]
}
greater_than_individual<- counts[rowSums(counts[,-1])>0, ] 
write.table(greater_than_individual,file="complete_virome_counts_unique_noTanRep_final.csv",sep=",")

