# Pseudospectral Integration Matrix (PSIM)

PSIM <- function(N,n,v,U1,U2) {

#---Chebyshev polynomials and their integrals---#
# Polynomials
TT <- cos( matrix(0:(N+2),N+3,1) %*% matrix(0:N,1,N+1) * pi/N )

# Integrals
dTm <- array(rep(0,(N+3)*(N+1)*(n+1)),c(N+3,N+1,n+1))
dTm[1:(N+3),1:(N+1),1] <- TT
for (m in 1:n) {
	dTm[1,1:(N+1),m+1] <- dTm[2,1:(N+1),m]
	dTm[2,1:(N+1),m+1] <- dTm[3,1:(N+1),m]/4
	for (k in 3:(N+3-m)) {
		dTm[k,1:(N+1),m+1] <- ( (dTm[k+1,1:(N+1),m]/k) - (dTm[k-1,1:(N+1),m]/(k-2)) )/2
	}
}
dTm <- aperm(dTm,c(2,1,3))

#---Beta coefficients---#
# Chebyshev weights
Cw      <- matrix(rep(1,N+1),1,N+1)
Cw[1]   <- 2
Cw[N+1] <- 2
Cw      <- N*(t(Cw) %*% Cw)

# Coefficients
b <- TT[1:(N+1),1:(N+1)] - TT[1:(N+1),v] %*% ( solve(TT[(N-n+2):(N+1),v],t(TT[(N-n+2):(N+1),1:(N+1)])) )
b <- 2*b/Cw

#---Boundary conditions---#
# Matrix system
Cm1 	   <- matrix(rep(0,n^2),n,n)
Cm1[1,1:n] <- rep(1,n)
if (n>=2) {
	for (k in 2:n) {
		Cm1[k,k:n] <- ( ( (k:n) - 1)^2 - (k-2)^2 ) / (2*k - 3)
		Cm1[k,k:n] <- Cm1[k,k:n]*Cm1[k-1,k:n]
	}
}
ind <- matrix(rep(0:(n-1),n),n,n)
Cm2 <- Cm1*(-1)^( t(ind) + ind )

# Rhs
dTm <- aperm(dTm,c(3,2,1))
P1  <- dTm[(n+1):2,1:(N+1),1]
P2  <- dTm[(n+1):2,1:(N+1),N+1]

# Free parameters
if (U1 == '') {
	RR <- solve( U2 %*% Cm2 )
	rr <- RR %*% ( U2 %*% P2 )
}
else if (U2 =='') {
	RR <- solve( U1 %*% Cm1 )
	rr <- RR %*% ( U1 %*% P1 )
}
else {
	RR <- solve( rbind( U1 %*% Cm1, U2 %*% Cm2 ) )
	rr <- RR %*% rbind( U1 %*% P1 , U2 %*% P2  )
}

#---Construction of PSIM---#
dTm <- aperm(dTm,c(3,2,1))
B   <- ( dTm[1:(N+1),1:(N+1),n+1] - dTm[1:(N+1),1:n,1] %*% rr ) %*% b
B[1:(N+1),v] <- dTm[1:(N+1),1:n,1] %*% RR
return(B)
}
