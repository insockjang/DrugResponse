require(RCytoscape)
require(graphite)

rCytoscapePlot<-function(drugID,pathwayName,curAlgorithm){
  
  
  cy = CytoscapeConnection()
  window.title = pathwayName #allPathways[k]
  if (window.title %in% as.character (getWindowList(cy)))
    deleteWindow (cy, window.title)
  cw = new.CytoscapeWindow (window.title, graphObject)
  
  #################### Node size change by xxx
  #############################################
  setNodeSizeDirect (cw, 'A', 32)
  # setNodeColorDirect (cw, 'A','#FFFF00' )
  displayGraph (cw)
  layoutNetwork (cw, 'jgraph-spring')
  redraw (cw)
  #################### Node size change by xxx
  #############################################
  
  #  send graph to Cytoscape.  this only needs to be done once!
  displayGraph (cw)
  
  #  ask Cytoscape to layout the graph
  layoutNetwork(cw, 'jgraph-spring')
  
  #  instruct Cytoscape to use each node's 'label' attribute as the value for the visible label it draws on the node
  # setNodeLabelRule (cw, 'label')
  #  experiment on your own, if you wish, with this occasionally useful rule
  # setNodeLabelRule (cw, 'nodeType')
  
  #  slighly more complicated rules.  phosphorylation is a directed activity, so put an arrow only on the substrate end
  #  binding is reciprocal, so put a circle on both ends of that edge
  
  setEdgeTargetArrowRule (cw, 'edgeType', c ('directed'), c ("Arrow"))
  setEdgeSourceArrowRule (cw, 'edgeType', c ('directed'), c ('None'))
  
  #  now ask Cytoscape to redraw the graph using these rules
  redraw (cw)
  
  #  set window and network sizes
  setWindowSize (cw, 1600, 1200)
  fitContent (cw)
  setZoom (cw, 0.8 * getZoom (cw))

  
}