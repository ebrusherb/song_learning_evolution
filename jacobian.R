####### song genetic, preference imprinted from mother

sigma2=1
sigmay2 = (5+2*sqrt(6))/3*sigma2+0
sigmax2 = Re(polyroot(c(2*sigma2^2,5*sigma2-3*sigmay2,3))[2])
J= matrix(c(1+sigmax2/4*(sigmay2-sigma2-4*sigmax2*sigmay2/(sigma2+sigmax2))/(sigma2+sigmax2)^2,sigmax2/2/(sigma2+sigmax2),1/2*sigma2*sigmay2/(sigma2+sigmax2)^2,1/2),nrow=2,byrow=TRUE)
print(eigen(J))