f = function(x,y,i,j){
	xy = matrix(c(x,y),2,1)
	ij = matrix(c(i,j),2,1)
	M  = matrix(c(5+i,2,2,1+j),2,2)
	return(0.5 * t(xy) %*% M %*% xy + t(xy) %*% ij - 6)
}
