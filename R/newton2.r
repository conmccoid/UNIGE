newton = function(x0,epsi,f,fp){
  err = abs(f(x0))
  while (err >= epsi && is.finite(err)){
    x0 = x0 - f(x0)/fp(x0)
    err = abs(f(x0))
  }
  return(x0)
}

f = function(x){return(x^3 - 4*x^2 - 1)}

fp= function(x){return(3*x^2 - 8*x)}

fpp = function(x){return(6*x - 8)}

print(newton(1,1e-8,fp,fpp))