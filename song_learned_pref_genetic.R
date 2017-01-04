sigma2 = 0.0
steps = 500

Rvals = seq(0,3,by=0.1)
Rvals = Rvals[-1]
xR = length(Rvals)
rhovals = seq(0,1,by=0.01)
xrho = length(rhovals)
sigmax2_init = 1

final_x = array(NA,dim=c(xR,xrho))
final_y = array(NA,dim=c(xR,xrho))
Fmat = array(NA,dim=c(xR,xrho))

for(i in 1:xR){
	for(j in 1:xrho){
		R = Rvals[i]
		rho = rhovals[j]
				Fmat[i,j] = sqrt(rho^2*(R^2-1)+2*(rho*R+1))
		sigmax2 = (sigmax2_init)^2
		sigmay2 = (R*sqrt(sigmax2_init))^2
		cov = rho*sqrt(sigmax2*sigmay2)
		for(t in 2:steps){
		rho = cov[t-1]/sqrt(sigmax2[t-1]*sigmay2[t-1])
		sigmax2[t] = sigmax2[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2
		sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2[t-1]))+1/4*cov[t-1]^2/(sigma2+sigmax2[t-1])^2*(sigmay2[t-1]-sigma2-sigmax2[t-1])
		cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2+1/2*sigmax2[t-1]*sigmay2[t-1]/(sigma2+sigmax2[t-1])
		}
		final_x[i,j]=sigmax2[steps]
		final_y[i,j]=sigmay2[steps]
	}
}

breaks=c(seq(min(final_x),1,length.out=5),seq(2,max(final_x),length.out=5))
layout(matrix(c(1,2),ncol=2),widths=c(3,1))
image(Rvals,rhovals,final_x,col=c(heat.colors(length(breaks)-1)),breaks=breaks);
image.scale(log(final_x),col=c(heat.colors(length(breaks)-1)),breaks=log(breaks),horiz=FALSE)

