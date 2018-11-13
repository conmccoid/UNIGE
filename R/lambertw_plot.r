# Plotting script for lambertw_fractal

filled.contour(x,y,log10(itersave),
	color.palette=rainbow,
	xlim=c(-l,l),ylim=c(-l,l),
	plot.axes = { points(Re(zsave),Im(zsave),col=2,pch=20) ; axis(1) ; axis(2) },
	main='Basin of attraction',xlab='Re(p)',ylab='Im(p)')

