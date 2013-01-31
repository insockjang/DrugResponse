require(Rgraphviz)
require(gplots)
require(plotrix)
pathwayLayout<-function(graphObject){
  # layout might be # “dot”, “neato”, “twopi”, “circo”, or “fdp”  
  g<-graphObject  
  
  # node shapes, sizes, and Fill color.
  Attrs <- list(graph=list(rankdir="LR"),node = list(fixedsize = TRUE),edge = list(arrowsize = 0.3))
  
  z1<-agopen(g,name="temp",layoutType = "fdp",attrs = Attrs)
  return(z1)
  
}

pathwayLayout1<-function(graphObject,features = NULL,corr){
  # layout might be # “dot”, “neato”, “twopi”, “circo”, or “fdp”  
  g<-graphObject  
  
  # selected feature 
  X<-sort(features,decreasing = T,index.return=TRUE)
  C<-corr[X$ix]
  Name<-names(X$x)
  Name1<-gsub("_expr","",Name)
  Name2<-gsub("_copy","",Name1)
  Name3<-gsub("_mut","",Name2)
  
  #   Colors<-colorRampPalette(brewer.pal(9,"RdBu"))(length(features))
  rw.palette <- colorRampPalette(c("red", "white"),space = "rgb")
  ColorsPos = rw.palette(length(features))
  bw.palette <- colorRampPalette(c("blue", "white"),space = "rgb")
  ColorsNeg = bw.palette(length(features))
  
  
  # Note: Fill, shape, width & height are all vectors which I set to get custom
  # node shapes, sizes, and Fill color.
  Attrs <- list(graph=list(rankdir="LR"))
  Attrs <- list(edge = list(arrowsize = 0.3))
  attrs <- list(node = list(fixedsize = TRUE))
  
  makeNodeDrawFunction <- function(x, y, nodeName) {
    force(x)
    force(y)
    force(nodeName)
    function(node, ur, attrs, radConv) {
      nc <- getNodeCenter(node)      
      names(x)<-colnames(x)      
      pieGlyph(x, xpos = getX(nc), ypos = getY(nc), radius = 25,col = y, cex = 1.5)
      text(getX(nc), getY(nc), nodeName, cex = 1.5, col = "black",font = 4)
    }
  }
  
  makeNodeDrawFunction1 <- function(x, y, nodeName) {
    force(x)
    force(y)
    force(nodeName)
    function(node, ur, attrs, radConv) {
      nc <- getNodeCenter(node)      
      names(x)<-colnames(x)      
      draw.circle(x=getX(nc),y=getY(nc),radius = 25,nv=100, border="black", col = y, lty=1, lwd=1)      
      text(getX(nc), getY(nc), nodeName, cex = 1.5, col = "black",font=4)
    }
  }
  
  drawFuns<-list()
  for(k in 1:numNodes(g)){
    
    nodeName<-nodes(g)[k]
    a<-which(is.na(match(Name3,nodeName))==0)
    
    if(length(a)<1){
      drawFuns[[k]]<-makeNodeDrawFunction1(1,0,nodeName)
    } else{
      name1<-c()
      cols<-c()
      for(kkk in 1:length(a)){
        if(length(grep("_expr",Name[a[kkk]]))!=0){name1<-c(name1,"E")}
        if(length(grep("_copy",Name[a[kkk]]))!=0){name1<-c(name1,"C")}
        if(length(grep("_mut",Name[a[kkk]]))!=0){name1<-c(name1,"M")}
        
        if(C[a[kkk]]>=0){
          cols<-c(cols,ColorsPos[a[kkk]])
        } else{
          cols<-c(cols,ColorsNeg[a[kkk]])
        }
      }
      counts<-matrix(1/length(a),ncol=length(a))
      colnames(counts)<-name1              
      drawFuns[[k]]<-makeNodeDrawFunction(counts,cols,nodeName)    
    }
  }
  
  return(drawFuns)
  
}

