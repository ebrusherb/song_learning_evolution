# sigma2 = 0.000
# steps = 50

# sigmax2 = c()
# sigmay2 = c()
# cov = c()

# sigmax2[1] = 0.01
# sigmay2[1] = 0.2
# rho = 0
# cov[1] = rho*sqrt(sigmax2*sigmay2)



# # # song learned from father, preference genetic:
# for(t in 2:steps){
	# rho = cov[t-1]/sqrt(sigmax2[t-1]*sigmay2[t-1])
	# sigmax2[t] = sigmax2[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2
	# sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2[t-1]))+1/4*cov[t-1]^2/(sigma2+sigmax2[t-1])^2*(sigmay2[t-1]-sigma2-sigmax2[t-1])
	# cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2+1/2*sigmax2[t-1]*sigmay2[t-1]/(sigma2+sigmax2[t-1])
# }

# par(mfrow=c(1,2));plot(sigmax2,ylim=c(0,3));points(sigmay2,col='red');plot(cov/sqrt(sigmax2*sigmay2))

# r = c(sqrt(sigmay2[1]/sigmax2[1]),exp(diff(log(sqrt(sigmay2)))))
# rho = cov/sqrt(sigmax2*sigmay2)
# P = matrix(0,nrow=steps)
# P[1] = rho[1]*sqrt(sigmay2[1]/sigmax2[1])
# P[2:steps] = 1/2^(1:(steps-1))*(P[1]-1)+1

# # #song learned obliquely, preference genetic
# # sigmay2 = c()
# # cov = c()

# # sigmax2 = 0.5
# # sigmay2[1] = 2
# # rho = 0.5
# # cov[1] = rho*sqrt(sigmax2*sigmay2)

# # for(t in 2:steps){
	# # sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2))+1/4*cov[t-1]^2/(sigma2+sigmax2)^2*(sigmay2[t-1]-sigma2-sigmax2)
	# # cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2)+sigmay2[t-1]*sigmax2)/(sigma2+sigmax2)^2+1/2*sigmax2*sigmay2[t-1]/(sigma2+sigmax2)
# # }

# # song genetic, preference genetic

# # sigma = 0.2899
# sigmay = 0.5

# poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# # sigmax = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# sigmax = poss 
# Q = (sigma^2*(sigma^2+sigmax^2)+sigmay^2*sigmax^2)/(sigma^2+sigmax^2)^2+2*sigmax*sigmay/(sigma^2+sigmax^2)+1
# print(sigmax)
# print(Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]))
# print(Q)
# v=seq(0,sigmay+0.2,by=0.01);plot(v,(sigma^2*(sigma^2+v^2)+sigmay^2*v^2)/(sigma^2+v^2)^2+2*v*sigmay/(sigma^2+v^2)+1,t='l');abline(h=4);abline(v=Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]),col='red')

par(mfrow=c(2,2))
# rho is already set at 1

sigma = 0.1
sigmay = 0.5
poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
solved = c(max(Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])),sigmay)
sigmax = solved[1]-0.01
Q = c()

for(t in 1:10){
	Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	sigmay = c(sigmay,Q[t]*sigmay[t]/4)
	sigmax = c(sigmax,Q[t]*sigmax[t]/4)
}

plot(sigmax);points(sigmay,col='red');abline(h=solved[1]);abline(h=solved[2],col='red')
plot(Q);abline(h=4)


sigma = 0.1
sigmay = 0.5
poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
solved = c(max(Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])),sigmay)
sigmax = solved[1]-0.01
Q = c()
rho = 0.8
cov = rho*sigmax*sigmay

for(t in 1:10){
	Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	sigmay = c(sigmay,Q[t]*sigmay[t]/4+(1-rho[t]^2)/4*sigmax)
	sigmax = c(sigmax,Q[t]*sigmax[t]/4)
}

plot(sigmax);points(sigmay,col='red');abline(h=solved[1]);abline(h=solved[2],col='red')
plot(Q);abline(h=4)

# rho evolves



# # sigma2 = 0.1
# rho = 1
# xvals = seq(0,5,by=0.05)
# yvals = seq(0,5,by=0.05)
# g =expand.grid(x=xvals,y=yvals)
# Q = with(g,sigma2*(sigma2+x)+x*y+2*rho*sqrt(x)*sqrt(y)*(sigma2+x)-3*(sigma2+x)^2)
# pdf('/Users/eleanorbrush/Documents/research/song_learning_evolution/Q_contour.pdf',width=3.4,height=3.4)
# contour(xvals,yvals,matrix(Q,nrow=length(xvals)),levels=seq(-1,1,length.out=7),xlab='Sigmax2',ylab='Sigmay2')
# abline(0,1,col='red')
# dev.off()