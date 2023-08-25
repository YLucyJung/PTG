# R v4.2.1
# chromatin interaction plot for CASR in parathyroid and other tissues
3d.interaction.plot<-function(){

  load("data/fig2.data.RData")
  # t: loops overlapping with promoters with a 5kb margin from hi-c parathyroid
  id=which(t$V11=="CASR"); # interactions with CARS promoter in PTG
  t[id,]

  # d: chromatin interactions in other tissues from pro-HiC, Jung et al.
  id2=which(d$Promoter=="CASR") # interacctions with CASR promoter in other tissues 
  d[id2,]

  #plot
  #ptg interaction first
  s=121760e3;e=122020e3;
  paste0("chr3:",s,"-",e)  
  TSS=121902530
  t3=data.frame(TSS,rowMeans(t[id,5:6]))
  
  layout(matrix(1:28),height=c(3,rep(1,times=27)),width=1,respect=F);  
  par(mar=c(2,4,2,10));
  plot(1,type="n",xlim=c(s,e), ylim=c(0,1), ylab=" ", xlab=" ", axes=F,xaxs = "i",yaxs = "i")
  abline(h=0)
  theta=seq(0,pi, len=100)
  
  for (i in seq(1,dim(t3)[1])){
  x = (abs(t3[i,1]-t3[i,2])/2)*cos(theta)+((t3[i,1]+t3[i,2])/2)
  y = sin(theta)*min((abs(t3[i,1]-t3[i,2])/2)/((e-s)/4),1)
  lines(x,y, xpd=TRUE,col="dodgerblue3")
  }
  mtext(3,text=paste0("chr3:",s,"-",e))
  mtext(4,text="Parathyroid",las=2)

  #draw interactions in other tissues  
  tissue=unique(d$Tissue_type)[1]
  tissues=unique(d$Tissue_type)
  for(tissue in tissues){
  print(tissue)
  id2=which(d$Tissue_type==tissue & d$Promoter=="CASR")
  par(mar=c(0.5,4,0,10),xpd=F);
  if(length(id2)==0){plot(1,type="n",xlim=c(s,e), ylim=c(0,1), ylab=" ", xlab=" ", axes=F,xaxs = "i",yaxs = "i");abline(h=0);mtext(4,text=tissue,las=2,cex=0.7)}
  if(length(id2)>0){
    t4=data.frame(TSS,unlist(lapply(strsplit(d[id2,2],"[.]"),function(dd) mean(as.numeric(dd[2:3])))))
  plot(1,type="n",xlim=c(s,e), ylim=c(0,1), ylab=" ", xlab=" ", axes=F,xaxs = "i",yaxs = "i")
                    
  abline(h=0);mtext(4,text=tissue,las=2,cex=0.7)
  theta=seq(0,pi, len=100)
  for (i in seq(1,dim(t3)[1])){
  x = (abs(t4[i,1]-t4[i,2])/2)*cos(theta)+((t4[i,1]+t4[i,2])/2)
  if(tissue=="Human ES cells" | tissue=="Right Ventricle"){y = sin(theta)*min((abs(t4[i,1]-t4[i,2])/2)/((e-s)/4),1)}
  else{y = sin(theta)*(abs(t4[i,1]-t4[i,2])/2)/((e-s)/4)}
  lines(x,y, xpd=F,col="dodgerblue3")
}
                  }  
}
}
