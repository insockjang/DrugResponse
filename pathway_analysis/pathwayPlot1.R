require(Rgraphviz)
require(gplots)
pathwayPlot1<-function(graphObject,features = NULL,outType="fdp",Fontsize=15,Width=10,Height=10,wGraph=100,hGraph=100){
  # layout might be # “dot”, “neato”, “twopi”, “circo”, or “fdp”  
  g<-graphObject  
  names(labels) = labels = nodes(g)
  #   fontscale = round(200/sqrt(length(labels)))  
  
  Fill <- sapply(labels, function(x) {0})
  
  if(!is.null(features)){
    # selected feature 
    Name<-names(features)
    X<-sort(features,index.return=TRUE)
    Name<-names(X$x)
    
    Name1<-strsplit(Name,"_")  
    Name2<-c()
    for(k in 1:length(Name)){
      Name2<-c(Name2,Name1[[k]][1])
    }  
    
    
#     Colors<-colorRampPalette(brewer.pal(9,"Greens"))(length(reference))
    Colors<-colorRampPalette(brewer.pal(9,"RdBu"))(length(reference))                                     
    namePathway<-names(labels)
    b<-match(namePathway,Name2)
    b[is.na(b)]<-0
  }
  else{
    b<-rep(0,length(Fill))
  }
  #   fixedsize = rep(FALSE, times=length(nodes(g)))
  
  for(k in 1:length(Fill)){
    if(b[k]==0){
      Fill[k] = "yellow"
    }
    else{
      Fill[k]=Colors[b[k]]
    }
  }
  
  Shape <- sapply(labels, function(x) {"circle"})  
   
  names(Fill) = names(Shape) = nodes(g)
  
  # Note: Fill, shape, width & height are all vectors which I set to get custom
  # node shapes, sizes, and Fill color.
  attrs <- list(graph=list(rankdir="LR"))
  attrs <- list(edge = list(arrowsize = 1.0))
  attrs <- list(node = list(fixedsize = TRUE))
  attrs <- list(graph = list(size = c(wGraph,hGraph)))
  x <- layoutGraph(g,layoutType = outType,attrs = attrs) 
  nodeRenderInfo(x) = list(shape = Shape, fill = Fill, height = Height, fontsize = Fontsize)
  renderGraph(x)                                         
}

