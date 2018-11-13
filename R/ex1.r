# Driver for PSIM example

source("psim.r")

# Grid
n <- 1
N <- 64
x <- cos( 0:N * pi / N )

# BCs
U2 <- ''
U1 <- 1

# PSIM
B <- PSIM(N,n,1,U1,U2)

# Solution
y <- B %*% matrix(rep(1,N+1),N+1,1)

plot(x,y,type="p",col=1)
lines(x,x,type="l",lty=2,col=2)
title(main="Test of PSIM",xlab="x",ylab="y")
legend(locator(1),legend=c("Numerical solution","Exact solution"),col=c(1,2),lty=1:2)
