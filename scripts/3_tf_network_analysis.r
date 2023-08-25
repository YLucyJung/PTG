
#R v.4
#enrichment comparison TSS-proximal vs TSS-distal for MAFB and GATA3
comparision.TSS.proximal.distall.mafb.gata3<-function(){
  tf=data.frame(proximal=c(913,1832),distal=c(1983,3629))
  #these numbers from predictied motif sites in TSS-poximal regions and TSS-distal regions
  rownames(tf)=c("MAFB","GATA3")
  b=data.frame(proximal=7531,distal=6842);
  (tf$proximal/rowSums(tf))/(b$proximal/rowSums(b)) #proximal
  #     MAFB     GATA3     
  #0.6016820 0.6402478  
  (tf$distal/rowSums(tf))/(b$distal/rowSums(b)) #distal
  #     MAFB     GATA3      
  #1.4384293 1.3959798    
  fisher.test(matrix(c(tf$proximal[1],tf$distal[1],b$proximal,b$distal),ncol=2))
  fisher.test(matrix(c(tf$proximal[2],tf$distal[2],b$proximal,b$distal),ncol=2))
  d=matrix(c(0.6016820,1.4384293,0.6402478,1.3959798),nrow=2)
  rownames(d)=c("Proximal","Distal");
  colnames(d)=c("MAFB","GATA3");
  t=barplot(d,ylim=c(0,1.6),beside=T,col=c("#fc9272","#9ecae1"),border=F,ylab="Enrichment");
  abline(h=1,lty="dashed")
  legend("topright",fill=c("#fc9272","#9ecae1"),legend=c("Proximal","Distal"),border=F,bty='n')
}

#TF network plot
tf.network.plot<-function(){
  # first we found TFs motifs at the GCM2 proximal regions (5kb window)
  # filter out silent genes in parathyroids or genes that are not associated with super-enhancers in parathyroids
  # For those genes, we made connections between genes and constructed network with nodes and links.
  # nodes and links for the tf network

  load("data/fig3.data.RData")
  head(nodes)
  head(links)
  nrow(nodes); length(unique(nodes$id))
  nrow(links); nrow(unique(links[,c("from", "to")]))

  library(igraph);
  net <- graph.data.frame(links, nodes, directed=T)
  x11();plot(net, edge.arrow.size=.4,vertex.frame.color="#ffffff");
    plot(net, layout=layout_in_circle,arrow.size=.2,vertex.frame.color="#ffffff");
  V(net)$color=c("red","blue");
    plot(net, layout=layout_in_circle,arrow.size=.2,vertex.frame.color="#ffffff");

  #ordered by connectivity
  genes=nodes$id;
  names(genes)=genes;
  counts=unlist(lapply(genes,function(gene){
    sum(links$from==gene)+sum(links$to==gene)}))
  order(counts)
  nodes2=genes[rev(order(counts))][c(2,1,3:15)]
  net2 <- graph.data.frame(links, nodes2, directed=T)
  
  library(RColorBrewer)
  cols=rev(unlist(lapply(1:9,function(i) {col=brewer.pal(n=9, name="Oranges")[i];
                                   return(c(col,col))})));
  V(net2)$color=cols[3:17]
  col2=rep("white",times=15);
  auto=unlist(lapply(nodes2,function(gene){
    id=which(links$from==gene & links$to==gene)
    if(length(id)>0){return(gene)}}))
  ids=match(auto,nodes2);
  col2[ids]="black"
  id=which(nodes2=="GATA3");
  col2[id]=brewer.pal(n=9, name="YlGn")[5];
  x11();
plot(net2, vertex.size=30,layout=layout_in_circle,arrow.size=.2,vertex.frame.color=col2, vertex.label.color=c("black"),vertex.label.family="Helvetica")
}
