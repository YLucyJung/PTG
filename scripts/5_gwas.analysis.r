# R v4.2.1
# Enrichment of GWAS SNPs of relevant phenotypes in regulatory elements of parathyroids


enrichment.gwas.snp.chromatin.state<-function(){
 load("data/f5.data.RData")
 
 #d:# of GWAS SNP in each regulatory regions such as chromatin states, DHS, super-enhancers, GCM2 peaks 
 
  #calculate p-value
  #compare with background
 
  #chromatin states
  #chromatin.states: chromatin state segmentations in PTG
  t=chromatin.states
  dist=lapply(1:18,function(i){
    id=which(t$V4==paste0("E",i));
    sum(t$V3[id]-t$V2[id])
  })
  unlist(dist)
  b.dist=unlist(dist); #background distribution 

  #p-value by fisher exact test
  p=list();
  for(j in 1:3){ #for each phenotypes
  o.dist=d[j,1:18]
  fisher.p=lapply(1:18,function(i){ #18 chromatin states
    fisher.test(matrix(c(o.dist[i],round(b.dist[i]/sum(b.dist)*sum(o.dist)),sum(o.dist),sum(o.dist)),nrow=2),alternative="greater")$p.value;
  })
  p1=unlist(fisher.p)[1:17]; 

  #DHS
  #dhs: DHS regions
  t=dhs;
  b=sum(t$V3-t$V2); 
  unique(t$V1) #it contains chrY

  #hg19:hg19 genome regions per chromosome
  t2=hg19
  g=sum(t2$V3[1:24])
  p2=fisher.test(matrix(c(d[j,19],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2),alternative="greater")$p.value;
  p2

  #se: super enhancer regions
  t=se
  b=sum(t$V3-t$V2); 
  #sort(unique(t$V1)) #it contains chrY
  p3=fisher.test(matrix(c(d[j,20],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2),alternative="greater")$p.value;
  p3

  #gcm2: GCM2 peaks
  t=gcm2
  b=sum(t$V3-t$V2); 
  sort(unique(t$V1)) #it contains chrM
  g=sum(t2$V3[1:25])
  p4=fisher.test(matrix(c(d[j,21],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2),alternative="greater")$p.value;
  p4  

  p[[j]]=c(p1,p2,p3,p4);
  p
  print(j);
  }
  po=do.call(rbind,p);

  p=p.adjust(po,method="BH") #adjust p-value
  p=matrix(p,nrow=3)
 
    
  #add odd ratio
  or=list();
  for(j in 1:3){
  #chromatin states
  o.dist=d[j,1:18]
  fisher.p=lapply(1:18,function(i){
    fisher.test(matrix(c(o.dist[i],round(b.dist[i]/sum(b.dist)*sum(o.dist)),sum(o.dist),sum(o.dist)),nrow=2))$estimate;
  })
  or1=unlist(fisher.p)[1:17];
  or1

  #DHS
  t=dhs;
  b=sum(t$V3-t$V2); 
  unique(t$V1) #it contains chrY
  t2=hg19;
  g=sum(t2$V3[1:24])
  or2=fisher.test(matrix(c(d[j,19],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2))$estimate;
  or2 

  #SE
  t=se
  b=sum(t$V3-t$V2); 
  #sort(unique(t$V1)) #it contains chrY
  g=sum(t2$V3[1:24])
  or3=fisher.test(matrix(c(d[j,20],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2))$estimate;
  or3

  #gcm2
  t=gcm2
  b=sum(t$V3-t$V2); 
  sort(unique(t$V1)) #it contains chrM
  g=sum(t2$V3[1:25])
  or4=fisher.test(matrix(c(d[j,21],round(b/g*d[j,22]),d[j,22],d[j,22]),nrow=2))$estimate;
  or4

  or[[j]]=c(or1,or2,or3,or4);

  print(j);
  }
  or=do.call(rbind,or);
  or

  range(-1*log10(p)) #0-24.22
  -1*log10(0.05) #1.30103
  or[or==Inf]=0
  or

  #figure
  #p-value for color
  #odd ratio for size
  require(ggplot2)
  
  pp=-1*log10(p);
  pp[pp>=4]=4
  colnames(pp)=colnames(d)[-c(18,22)]
  rownames(pp)=rownames(d);
  dataf=as.data.frame(as.table(pp))

  tmp=as.data.frame(as.table(or))
  tmp$Freq[tmp$Freq>=10]=10;
  names(dataf)=c("Type","RegulatoryElement","logP");
  dataf$OR=tmp$Freq;

  x11();
  ggplot(dataf, aes(y = factor(Type),
             x = factor(RegulatoryElement))) +       
  geom_point(aes(colour = logP, 
                   size =OR))  +  
   scale_color_gradient(low = "yellow",  
                       high = "red")+      
  scale_size(range = c(0, 15))+            
  theme_bw()
}                         

 
 
