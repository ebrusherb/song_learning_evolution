# sigma2 = 0.000
# steps = 50

# sigmax2 = c()
# sigmay2 = c()
# cov = c()

# c=1 
# C = 4
# rho = -1
# ratio= (-1+sqrt(-1+C+rho^2))/rho
# sigmax2[1] = 1*c
# sigmay2[1] = (ratio*sqrt(sigmax2))^2
# sigmay2[1] = 1*c

# cov[1] = rho*sqrt(sigmax2*sigmay2)


# # # song learned from father, preference genetic:
# for(t in 2:steps){
	# rho = cov[t-1]/sqrt(sigmax2[t-1]*sigmay2[t-1])
	# sigmax2[t] = sigmax2[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2
	# sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2[t-1]))+1/4*cov[t-1]^2/(sigma2+sigmax2[t-1])^2*(sigmay2[t-1]-sigma2-sigmax2[t-1])
	# cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2+1/2*sigmax2[t-1]*sigmay2[t-1]/(sigma2+sigmax2[t-1])
# }

# rho = cov/sqrt(sigmax2*sigmay2)
# R = sqrt(sigmay2/sigmax2)
# F = sqrt(rho^2*(R^2-1)+2*(rho*R+1))

# par(mfrow=c(1,3));
# plot(sigmax2,ylim=c(0,3));
# points(sigmay2,col='red');
# plot(rho)
# # plot(rho,R,t='o')
# plot(F)


# r = c(sqrt(sigmay2[1]/sigmax2[1]),exp(diff(log(sqrt(sigmay2)))))
# rho = cov/sqrt(sigmax2*sigmay2)
# P = matrix(0,nrow=steps)
# P[1] = rho[1]*sqrt(sigmay2[1]/sigmax2[1])
# P[2:steps] = 1/2^(1:(steps-1))*(P[1]-1)+1

# # # #song learned obliquely, preference genetic
# # # sigmay2 = c()
# # # cov = c()

# # # sigmax2 = 0.5
# # # sigmay2[1] = 2
# # # rho = 0.5
# # # cov[1] = rho*sqrt(sigmax2*sigmay2)

# # # for(t in 2:steps){
	# # # sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2))+1/4*cov[t-1]^2/(sigma2+sigmax2)^2*(sigmay2[t-1]-sigma2-sigmax2)
	# # # cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2)+sigmay2[t-1]*sigmax2)/(sigma2+sigmax2)^2+1/2*sigmax2*sigmay2[t-1]/(sigma2+sigmax2)
# # # }

# song genetic, preference genetic

# # # sigma = 0.2899
# # sigmay = 0.5

# # poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# # # sigmax = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# # sigmax = poss 
# # Q = (sigma^2*(sigma^2+sigmax^2)+sigmay^2*sigmax^2)/(sigma^2+sigmax^2)^2+2*sigmax*sigmay/(sigma^2+sigmax^2)+1
# # print(sigmax)
# # print(Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]))
# # print(Q)
# # v=seq(0,sigmay+0.2,by=0.01);plot(v,(sigma^2*(sigma^2+v^2)+sigmay^2*v^2)/(sigma^2+v^2)^2+2*v*sigmay/(sigma^2+v^2)+1,t='l');abline(h=4);abline(v=Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]),col='red')

# # # rho fixed at 1
# # sigma = 0.5
# # sigmay = 2
# # poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# # solved = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# # sigmax = solved[1]-0.01
# # Q = c()

# # for(t in 1:10){
	# # Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	# # sigmay = c(sigmay,sqrt(Q[t])/2*sigmay[t])
	# # sigmax = c(sigmax,sqrt(Q[t])/2*sigmax[t])
# # }
# # par(mfrow=c(1,2))
# # plot(sigmax);points(sigmay,col='red');abline(h=solved);abline(h=sigmay[1],col='red')
# # plot(Q);abline(h=4)

# rho = 1
# sigmaxvals = seq(0,5,by=0.05)
# sigmayvals = seq(0,5,by=0.05)
# g =expand.grid(x=sigmaxvals,y=sigmayvals)
# Q = with(g,(sigma^2*(sigma^2+x^2)+x^2*y^2)/(sigma^2+x^2)^2+2*rho*x*y/(sigma^2+x^2)+1)
# pdf('/Users/eleanorbrush/Documents/research/song_learning_evolution/Q_contour.pdf',width=3.4,height=3.4)
# contour(xvals,yvals,matrix(Q,nrow=length(xvals)),levels=seq(3,5,length.out=3),xlab='Sigmax2',ylab='Sigmay2')
# # abline(0,1,col='red')
# abline(h=sigmay[1],col='red')
# abline(v=solved,col='red')
# dev.off()

# rho evolves
sigma = 0.1
sigmay = 2
poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
solved = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
sigmax = solved[2]-0.00135
Q = c()
rho = 0.8
cov = rho*sigmax*sigmay

for(t in 1:20){
	Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	sigmay = c(sigmay,sqrt(Q[t]/4*sigmay[t]^2+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2*(sigma^2+sigmax[t]^2-sigmay[t]^2)/(sigma^2+sigmax[t]^2)^2))
	sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
	cov = c(cov,Q[t]/4*cov[t]+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2))
	rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
}
par(mfrow=c(1,3))
plot(sigmax);points(sigmay,col='red');abline(h=solved);abline(h=sigmay[1],col='red')
plot(rho)
plot(Q);abline(h=4)