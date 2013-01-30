require(Rgraphviz)
require(gplots)
require(plotrix)
pathwayPlot3<-function(graphObject,file,interestGenes = NULL,hubStart = NULL,hubEnd = NULL,...){
  # layout might be # “dot”, “neato”, “twopi”, “circo”, or “fdp”  
  g<-graphObject  
  
  # Note: Fill, shape, width & height are all vectors which I set to get custom
  # node shapes, sizes, and Fill color.
  Attrs <- list(graph=list(rankdir="LR"),node = list(fixedsize = TRUE),edge = list(arrowsize = 0.3))
  
  makeNodeDrawFunction <- function(x, y, nodeName) {
    force(x)
    force(y)
    force(nodeName)
    function(node, ur, attrs, radConv) {
      nc <- getNodeCenter(node)      
      names(x)<-colnames(x)
      #       pieGlyph(x, xpos = getX(nc), ypos = getY(nc), radius = getNodeRW(node),col = y, cex = 1.5)
      #       pieGlyph(x, xpos = getX(nc), ypos = getY(nc), radius = 25,col = y, cex = 1.5)
      draw.circle(x=getX(nc),y=getY(nc),radius = 20,nv=100,border="black", col=y, lty=1, lwd=1)      
      text(getX(nc), getY(nc), nodeName, cex = 1.5, col = "black", font = 4)
    }
  }
  
  
  drawFuns<-list()
  for(k in 1:numNodes(g)){    
    nodeName<-nodes(g)[k]
    if(is.na(match(nodeName,union(interestGenes,union(hubStart,hubEnd))))){
      drawFuns[[k]]<-makeNodeDrawFunction(1,"yellow",nodeName)    
    }else{
      if(!is.na(match(nodeName,hubStart))){
        drawFuns[[k]]<-makeNodeDrawFunction(1,"blue",nodeName)    
      }
      if(!is.na(match(nodeName,hubEnd))){
        drawFuns[[k]]<-makeNodeDrawFunction(1,"lightgreen",nodeName)    
      }
      if(!is.na(match(nodeName,setdiff(interestGenes,union(hubStart,hubEnd))))){
        drawFuns[[k]]<-makeNodeDrawFunction(1,"red",nodeName)    
      }  
    }
    
  }
  
  
  z1<-agopen(g,name="temp",layoutType = "fdp",attrs = Attrs)
  
  
  png(file,width= 1.5*z1@boundBox@upRight@x, height=1.5*z1@boundBox@upRight@y)
  plot(z1,drawNode = drawFuns)
  dev.off()
  
}
# need to implement set process (different node colors and edge colors etc)
