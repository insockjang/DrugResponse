source("~/COMPBIO/trunk/users/jang/pathway_analysis/preRankedTest.R")
require(ggplot2)
myPathwayAnalysis <- setRefClass(Class = "myPathwayAnalysis",                                 
                                  fields=c("gseaResult","fetResult"),
                                  methods = list(
                                    initialize = function(...){
                                      return(.self)
                                    },
                                    
                                    gsea = function(referenceSet,geneSet, np=1000, w=1){
                                      a1=paste(geneSet,"_copy",sep ="")
                                      a2=paste(geneSet,"_expr",sep ="")
                                      a3=paste(geneSet,"_mut",sep ="")
                                      geneSet <-union(a1,union(a2,a3))                                     
                                      geneSet <- intersect(geneSet, names(referenceSet))
                                      .self$gseaResult <- preRankedTest(referenceSet,geneSet, np = np, w = w)
                                    },
                                    
                                    fet = function(AllGenes,geneSet,testSet){ 
                                      #AllGenes = union(allGenes,union(geneSet,testSet))
                                      Mat2x2 <- mat.or.vec(2,2)
                                      Mat2x2[1,1] <- length(intersect(testSet,geneSet))
                                      Mat2x2[2,1] <- length(setdiff(testSet,geneSet))
                                      Mat2x2[1,2] <- length(setdiff(geneSet,testSet))
                                      Mat2x2[2,2] <- length(AllGenes) - Mat2x2[1,1]- Mat2x2[1,2]- Mat2x2[2,1]
                                      
                                      .self$fetResult<-fisher.test(Mat2x2)                               
                                    },
                                    gseaPlot = function(reference,geneSet){
                                      a1=paste(geneSet,"_copy",sep ="")
                                      a2=paste(geneSet,"_expr",sep ="")
                                      a3=paste(geneSet,"_mut",sep ="")
                                      geneSet <- union(a1,union(a2,a3))                                     
                                      geneset <- intersect(geneSet, names(reference))                                      
                                      
                                      pn      <- matrix(1,ncol = length(reference));
                                      
                                      sortRef <- sort(reference, decreasing = TRUE)
                                      
                                      isInGeneset <- sign(match(names(sortRef), geneset, nomatch=0))
                                      
                                      # compute ES
                                      score_hit   <- cumsum((abs(sortRef * isInGeneset)))
                                      
                                      if (sum(score_hit) == 0){
                                        print("no significant genes in geneset. Returning ES = 0")
                                        break
                                      }
                                      
                                      score_hit   <- score_hit/score_hit[length(score_hit)]
                                      score_miss  <- cumsum( 1 - isInGeneset)
                                      score_miss  <- score_miss/score_miss[length(score_miss)]
                                      es_all      <- score_hit - score_miss
                                      
                                      ############################### Enrichment plot
                                      N<-length(reference)
                                      M<-length(geneset)
                                      ind<-1:N                                     
                                      
                                      pos.max.ES<-which.max(es_all)
                                      
                                      min.RES <- min(es_all)
                                      max.RES <- max(es_all)
                                      
                                      if (max.RES < 0.3) max.RES <- 0.3
                                      if (min.RES > -0.3) min.RES <- -0.3
                                      
                                      delta <- (max.RES - min.RES)*0.50
                                      min.plot <- min.RES - 2*delta
                                      max.plot <- max.RES
                                      max.corr <- max(reference)
                                      min.corr <- min(reference)
                                      Obs.correl.vector.norm <- (reference - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
                                      zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
                                      col <- ifelse(.self$gseaResult$es > 0, 2, 4)
                                      
                                      # Running enrichment plot
                                      
                                      # sub.string <- paste("Number of genes: ", N, " (in list), ", M, " (in gene set)", sep = "", collapse="")  
                                      # main.string <- paste("Gene Set :", PathwayName)
                                      sub.string1 <- paste("Number of genes: ", N, " (in list), ", M, " (in gene set)",sep = "", collapse="")  
                                      sub.string2 <- paste("Statistics of gene set: ", sprintf("%.3f",.self$gseaResult$p.value), " (normial Pvalue), ", sprintf("%.3f",.self$gseaResult$nes), " (normalized enrichment score)",sep = "", collapse="")  
                                      
                                      plot(ind, es_all, sub = sub.string2, xlab = sub.string1, ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = .5, col = col)
                                      
                                      
                                      for (j in seq(1, N)) {
                                        lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
                                      }
                                      
                                      lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
                                      lines(c(pos.max.ES, pos.max.ES), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
                                      for (j in 1:N) {
                                        if (isInGeneset[j] == 1) {
                                          lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
                                        }
                                      }
                                      
                                      lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
                                      lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
                                      
                                      
                                    }
                                    
                                    
                                    
                                    )
                                  
                                  )
