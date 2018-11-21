newton = function(x0,epsi,f,fp){
error = 1
while (error > epsi){
	x1 = x0 - f(x0)/fp(x0)
	error = abs(x1 - x0)
	x0 = x1
	}
return(x0)
}

f = function(x){
return(x^3 - 4*x^2 -1)
}

fp = function(x){
return(3*x^2 - 8*x)
}
