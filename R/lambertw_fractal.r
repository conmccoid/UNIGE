# The object of this file is to form the fractal pattern generated by Newton's method solving
# the nonlinear equation x e^x - 1 = 0 in the complex plane.

# Form the complex plane
l <- 5
d <- 2*l/500
x <- seq(-l,l,d)
y <- seq(-l,l,d)
N <- length(x)
M <- length(y)

# Initialize the variables
itersave <- matrix(rep(0,N*M),N,M)
zsave	 <- itersave
itermax  <- 1000
tol	 <- 1e-8

# Newton's iteration
for (i1 in 1:N) {
	for (i2 in 1:M) {
		z    <- complex(real = x[i1], imaginary = y[i2])
		error<- 1
		iter <- 0
		while (is.finite(error) && error > tol && iter < itermax) {
			iter <- iter + 1
			z_new<- (z^2 + exp(-z))/(z+1)
			error<- Mod(z_new-z)
			z    <- z_new
		}
		
		if (iter < itermax & is.finite(z)) {
			itersave[i1,i2] <- iter
			zsave[i1,i2]    <- z
		}
		else {
			itersave[i1,i2] <- NaN
			zsave[i1,i2]    <- NaN
		}
	}
}

# Plot contour map of number of iterations
filled.contour(x,y,log10(itersave),
	color.palette=rainbow,
	xlim=c(-l,l),ylim=c(-l,l),
	plot.axes = { points(Re(zsave),Im(zsave),col=2,pch=20) ; axis(1) ; axis(2) },
	main='Basin of attraction',xlab='Re(p)',ylab='Im(p)')
