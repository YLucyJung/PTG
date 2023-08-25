
# This code was run R v4.2.1

# tsne plot for ptg expression along with other tissues from GTEX
tsne.tissue.expression.plot<-function(){
 library(Rtsne)
 load("data/f1.data.RData")
 # a: data from GTEx_v7_Annotations_SampleAttributesDS.txt, containing SAMPID and SMTS
 types=unique(a$SMTSD);
 length(types) # number of tissue types: 54

 # e: GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct
 # containing expression levels in TPM for tissues and samples 

 zscore=function(n=n){
  z=(n-mean(n))/sd(n);
  return(z);
 }

  # add median
  # d: expression data from table cotaining gene expression levels in TPM from parathyroids, 8 individuals
  medians=apply(d[,3:10],1,median);
  d$median=medians;

  #t: medians of TPM per gene per GTEX tissues
  genes=t$genes;
  id=match(genes,d$Symbol);
  t$Parathyroid=d$median[id];
  z=apply(t[,2:55],1,zscore);
  zz=z[54,]
  zz=list();
  zz$genes=t$genes;
  zz$interZ=z[54,];
  zz$genes[zz$interZ>7] #highly parathryoid specific genes...GCM2, CASR, PTH included...and PAX1
  z=zscore(t[,55])
  zz$intraZ=z;
  zz$genes[zz$interZ>7 & zz$intraZ>0];
  genes1=zz$genes[zz$interZ>7 & zz$intraZ>0]; #parathyroid specific genes

 exe=list();
 tissues=setdiff(types,"Cells - Leukemia cell line (CML)")
 tissues=c(tissues,"Parathyroid")
 for(gene in genes1){
   print(gene);
  for(tissue in tissues){
  if(tissue=="Parathyroid"){
    id=which(d$Symbol==gene);
    exe[[gene]][[tissue]]=as.numeric(d[id,3:10]);
  }
  if(tissue != "Parathyroid"){
  id=which(a$SMTSD==tissue);
  samid=a$SAMPID[id];
  samid=gsub("-",".",samid);
  id2=match(samid,names(e));
  id2=id2[!is.na(id2)]
    id3=which(e$Description==gene)[1];
    exe[[gene]][[tissue]]=as.numeric(e[id3,id2]);
  }}}
  exe2=lapply(exe,function(d){
    unlist(d)})
  exe3=do.call(rbind,exe2)
  tsne_results=Rtsne(t(exe3),perplexity=60,check_duplicates=F);
  #this result can be slightly changed each time due to randomness

  library(randomcoloR)
  cols <- randomcoloR::distinctColorPalette(k = 54) #this result can be slightly changed each time due to randomness
  col.list=unlist(lapply(as.factor(gsub('[0-9]+', '',colnames(exe3))),function(d){
  cols[which(tissues==d)]}))
  tissue.type=gsub('[0-9]+', '',colnames(exe3))

  x11();
  par(mfrow=c(1,3))
  plot(tsne_results$Y,col="white",bg=col.list, pch = 21, cex = 1.4,xlab="tSNE1",ylab="tSNE2",axes=F,useRaster=T); #other tissues
  axis(1);axis(2,las=2)
points(jitter(tsne_results$Y[11689:11696,],amount=1),col = "black",bg="orange",pch=21,cex=1.6) #parathyroid tissues, added some jitter for visualization 
  plot.new()
  legend("topleft",col=cols[1:27],pch=19,cex=1.2,legend=tissues[1:27],bty='n')
  plot.new()
  legend("topleft",col=c(cols[28:53],"orange"),pch=19,cex=1.2,legend=tissues[28:54],bty='n')

  #tissue type legend
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '',xlim=range(tsne_results$Y),ylim=range(tsne_results$Y),xlab="tSNE1",ylab="tSNE2",axes=F);
  axis(1);
  axis(2,las=2)
  plot.new();
  plot.new()
}
