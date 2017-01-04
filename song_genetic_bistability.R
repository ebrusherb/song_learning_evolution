sigma = 0.5
rho_init=0
steps = 25

sigmaxvals = seq(0,5,length.out=30)
xx = length(sigmaxvals)
sigmayvals = seq(0,5,length.out=30)
xy = length(sigmayvals)
rhovals = seq(0,1,length.out=20)
xr = length(rhovals)

final_x_prefgen = array(NA,dim=c(xx,xy,xr))
final_y_prefgen = array(NA,dim=c(xx,xy,xr))
final_rho_prefgen = array(NA,dim=c(xx,xy,xr))

final_x_preffath = array(NA,dim=c(xx,xy,xr))
final_y_preffath = array(NA,dim=c(xx,xy,xr))
final_rho_preffath = array(NA,dim=c(xx,xy,xr))

final_x_prefmoth = array(NA,dim=c(xx,xy,xr))
final_y_prefmoth = array(NA,dim=c(xx,xy,xr))
final_rho_prefmoth = array(NA,dim=c(xx,xy,xr))

for(i in 1:xx){
	for(j in 1:xy){
		for(k in 1:xr){
		sigmax_init = sigmaxvals[i]
		sigmay_init = sigmayvals[j]
		rho_init = rhovals[k]
		
		sigmax = sigmax_init
		sigmay = sigmay_init
		rho = rho_init
		cov = rho_init*sigmax*sigmay
		Q = c()
		for(t in 1:steps){
			Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
		sigmay = c(sigmay,sqrt(Q[t]/4*sigmay[t]^2+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2*(sigma^2+sigmax[t]^2-sigmay[t]^2)/(sigma^2+sigmax[t]^2)^2))
		sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
		cov = c(cov,Q[t]/4*cov[t]+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2))
		rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
		}
		final_x_prefgen[i,j,k]=sigmax[steps]
		final_y_prefgen[i,j,k]=sigmay[steps]
		final_rho_prefgen[i,j,k]=rho[steps]
		
		sigmax = sigmax_init
		sigmay = sigmay_init
		rho = rho_init
		cov = rho_init*sigmax*sigmay
		Q = c()
		for(t in 1:steps){
			if(sigma!=0){
			Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
			sigmay = c(sigmay,sqrt(sigmax[t]^2*(sigma^2/(sigma^2+sigmax[t]^2)+sigmax[t]^2*sigmay[t]^2/(sigma2+sigmax[t]^2)^2)))
			sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
			cov = c(cov,1/2*sigmax[t]^2*(sigma^2/(sigma^2+sigmax[t]^2)+sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2)^2+cov[t]))
			rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
		} 
		else {
			Q = c(Q,(sigmay[t]^2)/(sigmax[t]^2)+2*rho[t]*sigmay[t]/(sigmax[t])+1)
			sigmay = c(sigmay,sqrt(sigmax[t]^2*(sigmay[t]^2/(sigmax[t]^2))))
			sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
			cov = c(cov,1/2*sigmax[t]^2*((sigmay[t]^2)/(sigmax[t]^2)+rho[t]*sigmay[t]/(sigmax[t])))
			rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
			}
		}
		final_x_preffath[i,j,k]=sigmax[steps]
		final_y_preffath[i,j,k]=sigmay[steps]
		final_rho_preffath[i,j,k]=rho[steps]
		
		sigmax = sigmax_init
		sigmay = sigmay_init
		rho = rho_init
		cov = rho_init*sigmax*sigmay
		Q = c()
		for(t in 1:steps){			
			Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
		sigmay = c(sigmay,sigmay_init)
		sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
		cov = c(cov,1/2*sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2)+1/2*cov[t])
		rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])		
		}
		final_x_prefmoth[i,j,k]=sigmax[steps]
		final_y_prefmoth[i,j,k]=sigmay[steps]
		final_rho_prefmoth[i,j,k]=rho[steps]
		}
	}
}

xbreaks=c(0,1e-10,seq(.1,10,length.out=20),1000)
ybreaks=c(0,1e-10,seq(.1,10,length.out=20),1000)
rhobreaks=seq(0,1,length.out=50)

layout(matrix(1:6,ncol=3))
k=10
j=9
image(sigmaxvals,sigmayvals,final_x_prefgen[,,k],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))
image(sigmaxvals,rhovals,final_x_prefgen[,j,],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))
image(sigmaxvals,sigmayvals,final_x_preffath[,,k],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))
image(sigmaxvals,rhovals,final_x_preffath[,j,],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))
image(sigmaxvals,sigmayvals,final_x_prefmoth[,,k],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))
image(sigmaxvals,rhovals,final_x_prefmoth[,j,],breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))