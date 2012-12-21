which2d <-function(mat,sym = TRUE){
	nRow <-dim(mat)[1]
	nCol <-dim(mat)[2]
	MAT <- c()
	if (sym == TRUE){
		for (i in 1:nRow){
			for (j in i:nCol){
				if( mat[i,j] == 1){
					MAT<-rbind(MAT,c(i,j))
				}			
			}
		}
	}
	else{
		for (i in 1:nRow){
			for (j in 1:nCol){
				if( mat[i,j] == 1){
					MAT<-rbind(MAT,c(i,j))
				}			
			}
		}
	}
	return(MAT)
}