# R v4.2.1
# chromatin environments(chromatin states) and expression states of PTG-specific genes in comparison of other tissues

ptg.gene.chromatin.expression.state<-function(){

  cols=c("255,0,0","255,69,0","255,69,0","255,69,0","0,128,0","0,100,0","194,225,5","194,225,5","255,195,77","255,195,77","255,255,0","102,205,170","138,145,208","205,92,92","189,183,107","128,128,128","192,192,192","255,255,255");
  #colors for each state. color codes are the same as in Roadmap Epigenomics

  cols=unlist(lapply(cols,function(d){
  dd=strsplit(d,",")[[1]];
  d1=as.numeric(dd[1]);
  d2=as.numeric(dd[2]);
  d3=as.numeric(dd[3]);
  rgb(d1,d2,d3,max=255)
})) #color changes into rbg

#make matrix for the gene
#each bin for 200bp
  load("data/f4.data.RData")
chromatin.matrix<-function(chr=chr,s=s,e=e,gene=gene,cols=cols,d=d,t=t){
  #heatmap for parathyroid
  #t: 18 chromatin states for parathyroid
if(sum(t$V1==chr & t$V2==s)>0){
  id1=which(t$V1==chr & t$V2==s)} else {id1=min(which(t$V1==chr & t$V2>=s))-1};
if(sum(t$V1==chr & t$V3==e)>0) {
  id2=which(t$V1==chr & t$V3==e);} else (id2=max(which(t$V1==chr & t$V3<=e))+1);
t.s=t[id1:id2,];

#expand table with 200bp
t2=lapply(1:dim(t.s)[1],function(i){
  dd=t.s[i,];
  if(i==1){d3=data.frame(dd$V1,seq(s,dd$V3-200,by=200),seq(s+200,dd$V3,by=200),dd$V4)}
  if(i==dim(t.s)[1]){d3=data.frame(dd$V1,seq(dd$V2,e-200,by=200),seq(dd$V2+200,e,by=200),dd$V4)}  
  if(!(i==1 | i==dim(t.s)[1])){d3=data.frame(dd$V1,seq(dd$V2,dd$V3-200,by=200),seq(dd$V2+200,dd$V3,by=200),dd$V4)}
  names(d3)=names(t.s);
  return(d3)
})
t3=do.call(rbind,t2);
names(t3)=names(t)
dim(t3)
states=as.numeric(gsub("E","",t3$V4));
mat=list();
mat[["PTG"]]=states; #chromatin matrix for parathyroid

#chromatin states for other tissues
#d: all 18 chromatin states for other tissues from Roadmap Epigenomics
tissues=names(d)  
tissues.s=tissues;

for(tissue in tissues.s){  
t=d[[tissue]]
print(tissue);

#make subset table
if(sum(t$V1==chr & t$V2==s)>0){
  id1=which(t$V1==chr & t$V2==s)} else {id1=min(which(t$V1==chr & t$V2>=s))-1};
if(sum(t$V1==chr & t$V3==e)>0) {
  id2=which(t$V1==chr & t$V3==e);} else (id2=max(which(t$V1==chr & t$V3<=e))+1);
t.s=t[id1:id2,];

#expand table with 200bp
t2=lapply(1:dim(t.s)[1],function(i){
  dd=t.s[i,];
  if(i==1){d3=data.frame(dd$V1,seq(s,dd$V3-200,by=200),seq(s+200,dd$V3,by=200),dd$V4,dd$V5,dd$V6,dd$V7,dd$V8,dd$V9)}
  if(i==dim(t.s)[1]){d3=data.frame(dd$V1,seq(dd$V2,e-200,by=200),seq(dd$V2+200,e,by=200),dd$V4,dd$V5,dd$V6,dd$V7,dd$V8,dd$V9)}  
  if(!(i==1 | i==dim(t.s)[1])){d3=data.frame(dd$V1,seq(dd$V2,dd$V3-200,by=200),seq(dd$V2+200,dd$V3,by=200),dd$V4,dd$V5,dd$V6,dd$V7,dd$V8,dd$V9)}
  names(d3)=names(t.s);
  return(d3)
})
t3=do.call(rbind,t2);
names(t3)=names(t)
dim(t3)

test=as.character(t3$V4);
states=unlist(lapply(test,function(d){
  as.numeric(strsplit(d,"_")[[1]][1])}))
mat[[tissue]]=states;
}

n=length(mat[[1]]); 
mat.sub=lapply(mat,function(d) if(length(d)==n)d);

mat.com=do.call(rbind,mat.sub);
  
return(mat.com)
}

#make plots
library(Sushi)

#Sushi_benes.bed: gene info in hg19
  
gene="DNAH11" #example 
x11();  
layout(matrix(1:3,nrow=3),height=c(0.7,0.5,3),width=1,respect=F);
par(mar=c(0.3,4,3,4));
gene=gene;
id4=which(Sushi_genes.bed$gene==gene);
id4
bu.genes=Sushi_genes.bed;
t7=bu.genes[id4,];
t7

    info=list();
    print(gene); 
    #g.t: gene table from hg19.refseq
    id=which(g.t$V4==gene);
    chr=g.t$V1[id][1];g.s=min(g.t$V2[id]);g.e=max(g.t$V3[id]);
    d=g.e-g.s;s=g.s-d*2/3;e=g.e+d*2/3;s=floor(s/200)*200;e=floor(e/200)*200;
    info[[gene]]=data.frame(gene=gene,s=s,e=e);#}
 

s=info[[gene]]$s;
e=info[[gene]]$e;
pg = plotGenes(t7,chrom=chr,chromstart=s,chromend=e ,packrow=F,maxrows=0.3,bheight=0.3,types="exon",plotgenetype="box",bentline=FALSE,labeloffset=.4,fontsize=0.8,arrowlength = 0.025,labeltext=TRUE)
mtext(3,at=s+(e-s)*0.1,text=paste0(chr,":",s,"-",e),cex=0.8);
box();
par(mar=c(0.3,4,0.3,4));
  
#panel for PTG
load("data/f4.data.RData")
mat.com=chromatin.matrix(chr=chr,s=s,e=e,gene=gene,cols=cols,d=d,t=t)
image((as.matrix(mat.com[1,])),col=cols,xaxs="i",yaxs="i",axes=F,useRaster=T);box()
mtext(side=4,text="Parathyroid",line=0.7)

#panel for other tissues
par(mar=c(3,4,0.3,4));
image(t(as.matrix(mat.com[-1,])),col=cols,xaxs="i",yaxs="i",axes=F,useRaster=T);box()
mtext(side=4,text="Other tissues",line=0.7)  
}
  

