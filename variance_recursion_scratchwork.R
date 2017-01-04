source('/Users/eleanorbrush/Desktop/image.scale.R')
### # # # song learned from father, preference genetic:
# sigma2 = 0.1
# steps = 50

# # c=1
# C = 4
# rho = 0.9
# ratio= (-1+sqrt(-1+C+rho^2))/rho
# sigmax2 = 1*c
# sigmay2 = (ratio*sqrt(sigmax2))^2
# sigmay2 = 1.1*c

# cov = rho*sqrt(sigmax2*sigmay2)
# for(t in 2:steps){
	# rho = cov[t-1]/sqrt(sigmax2[t-1]*sigmay2[t-1])
	# sigmax2[t] = sigmax2[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2
	# sigmay2[t] = 1/2*sigmay2[t-1]*(1+cov[t-1]/(sigma2+sigmax2[t-1]))+1/4*cov[t-1]^2/(sigma2+sigmax2[t-1])^2*(sigmay2[t-1]-sigma2-sigmax2[t-1])
	# cov[t] = 1/2*cov[t-1]*(sigma2*(sigma2+sigmax2[t-1])+sigmay2[t-1]*sigmax2[t-1])/(sigma2+sigmax2[t-1])^2+1/2*sigmax2[t-1]*sigmay2[t-1]/(sigma2+sigmax2[t-1])
# }

# rho = cov/sqrt(sigmax2*sigmay2)
# R = sqrt(sigmay2/sigmax2)
# F = sqrt(rho^2*(R^2-1)+2*(rho*R+1))

# # layout(matrix(1:4,ncol=2))
# plot(log(sigmax2),ylim=c(-3,1.5))#,ylim=c(0,3));
# points(log(sigmay2),col='red');
# plot(rho,ylim=c(0,1))
# # plot(rho,R,t='o')
# plot(F,ylim=c(0,2))
# plot(R)

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

##### # song genetic, preference genetic

sigma = 0.2899
sigmay = 0.5

# ##### finding sigmax that makes Q=4
# poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# # sigmax = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# sigmax = poss 
# Q = (sigma^2*(sigma^2+sigmax^2)+sigmay^2*sigmax^2)/(sigma^2+sigmax^2)^2+2*sigmax*sigmay/(sigma^2+sigmax^2)+1
# # print(sigmax)
# # print(Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]))
# # print(Q)
# v=seq(0,sigmay+0.2,by=0.01);plot(v,(sigma^2*(sigma^2+v^2)+sigmay^2*v^2)/(sigma^2+v^2)^2+2*v*sigmay/(sigma^2+v^2)+1,t='l');abline(h=4);abline(v=Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))]),col='red')

# #### rho fixed at 1, evolve sigmax,sigmay: sigmax such that Q=4 is an unstable equilibrium
# sigma = 0.5
# sigmay = 2
# poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# solved = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# sigmax = solved[2]-.001
# Q = c()

# for(t in 1:10){
	# Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	# sigmay = c(sigmay,sqrt(Q[t])/2*sigmay[t])
	# sigmax = c(sigmax,sqrt(Q[t])/2*sigmax[t])
# }
# par(mfrow=c(1,2))
# plot(sigmax);points(sigmay,col='red');
# abline(h=solved);
# abline(h=sigmay[1],col='red')
# plot(Q);abline(h=4)

# #### contour of Q as a function of sigmax,sigmay
# source('/Users/eleanorbrush/Desktop/image.scale.R')
# sigma = 0.5
# rho = 0.5
# sigmaxvals = seq(0,5,by=0.05)
# sigmayvals = seq(0,5,by=0.05)
# g =expand.grid(x=sigmaxvals,y=sigmayvals)
# Q = with(g,(sigma^2*(sigma^2+x^2)+x^2*y^2)/(sigma^2+x^2)^2+2*rho*x*y/(sigma^2+x^2)+1)
# # # # # # pdf('/Users/eleanorbrush/Documents/research/song_learning_evolution/Q_contour.pdf',width=3.4,height=3.4)
# # contour(sigmaxvals,sigmayvals,matrix(Q,nrow=length(sigmaxvals)),levels=seq(3,5,length.out=3),xlab='Sigmax2',ylab='Sigmay2')
# # # abline(0,1,col='red')
# # abline(h=sigmay[1],col='red')
# # abline(v=solved,col='red')
# # abline(0,1,col='blue')
# layout(matrix(c(1,2),ncol=2))
# breaks=c(0,4,seq(5,100,length.out=5))
# image(sigmaxvals,sigmayvals,matrix(Q,nrow=length(sigmaxvals)),breaks=breaks,col=c('black',heat.colors(length(breaks)-2)));
# abline(0,1,col='blue');
# abline(1,1,col='blue')
# image.scale(Q,breaks=breaks,horiz=FALSE,col=c('black',heat.colors(length(breaks)-2)))
# # # # # # # dev.off()

# # # # rho evolves
# sigma = 0.5
# sigmay = 1.5
# poss = polyroot(c(-2*sigma^4,2*sigmay*sigma^2,-5*sigma^2+sigmay^2,2*sigmay,-3))
# solved = Re(poss[intersect(which(Re(poss)>0),which(abs(Im(poss))<1e-14))])
# sigmax = solved[2]-0.00135
# sigmax = solved[2]-0.02905#116#29995
# sigmax = 0.5
# Q = c()
# rho = 0
# cov = rho*sigmax*sigmay

# for(t in 1:500){
	# Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
	# sigmay = c(sigmay,sqrt(Q[t]/4*sigmay[t]^2+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2*(sigma^2+sigmax[t]^2-sigmay[t]^2)/(sigma^2+sigmax[t]^2)^2))
	# sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
	# cov = c(cov,Q[t]/4*cov[t]+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2))
	# rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
# }
# par(mfrow=c(1,3))
# plot(sigmax);points(sigmay,col='red');abline(h=solved);abline(h=sigmay[1],col='red')
# plot(rho)
# plot(Q);abline(h=4)

##### this chunk can all be run together: 
# # look at bistability 
sigma = 0.5
rho_init=0.7
steps = 200

sigmaxvals = seq(0,10,length.out=50)
xx = length(sigmaxvals)
sigmayvals = seq(0,10,length.out=50)
xy = length(sigmayvals)

final_x = array(NA,dim=c(xx,xy))
final_y = array(NA,dim=c(xx,xy))
final_rho = array(NA,dim=c(xx,xy))

for(i in 1:xx){
	for(j in 1:xy){
		sigmax = sigmaxvals[i]
		sigmay = sigmayvals[j]
		Q = c()
		rho = rho_init
		cov = rho_init*sigmax*sigmay
		for(t in 1:steps){
			Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
		sigmay = c(sigmay,sqrt(Q[t]/4*sigmay[t]^2+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2*(sigma^2+sigmax[t]^2-sigmay[t]^2)/(sigma^2+sigmax[t]^2)^2))
		sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
		cov = c(cov,Q[t]/4*cov[t]+(1-rho[t]^2)/4*sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2))
		rho = c(rho,cov[t+1]/sigmax[t+1]/sigmay[t+1])
		}
		final_x[i,j]=sigmax[steps]
		final_y[i,j]=sigmay[steps]
		final_rho[i,j]=rho[steps]
	}
}

layout(matrix(c(1:6),ncol=2,byrow=TRUE));
xbreaks=c(0,seq(0.01,100,length.out=20))
image(sigmaxvals,sigmayvals,final_x,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)));
image.scale(final_x,horiz=FALSE,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))

ybreaks=c(0,seq(0.1,10,length.out=20))
image(sigmaxvals,sigmayvals,final_y,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)));
image.scale(final_y,horiz=FALSE,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)))

rhobreaks=c(0,seq(0.5,0.99,length.out=20),1)
image(sigmaxvals,sigmayvals,final_rho,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)));
image.scale(final_rho,horiz=FALSE,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)))

######  song genetic, preference from father

# sigma = 0.25
# rho_init=0
# steps = 50

# sigmaxvals = seq(0,10,length.out=50)
# xx = length(sigmaxvals)
# sigmayvals = seq(0,10,length.out=50)
# xy = length(sigmayvals)

final_x_v2 = array(NA,dim=c(xx,xy))
final_y_v2 = array(NA,dim=c(xx,xy))
final_rho_v2 = array(NA,dim=c(xx,xy))

for(i in 1:xx){
	for(j in 1:xy){
		sigmax = sigmaxvals[i]
		sigmay = sigmayvals[j]
		Q = c()
		rho = rho_init
		cov = rho_init*sigmax*sigmay
		for(t in 1:steps){
			if(sigma!=0){
			Q = c(Q,(sigma^2*(sigma^2+sigmax[t]^2)+sigmay[t]^2*sigmax[t]^2)/(sigma^2+sigmax[t]^2)^2+2*rho[t]*sigmax[t]*sigmay[t]/(sigma^2+sigmax[t]^2)+1)
		sigmay = c(sigmay,sqrt(sigmax[t]^2*(sigma^2/(sigma^2+sigmax[t]^2)+sigmax[t]^2*sigmay[t]^2/(sigma2+sigmax[t]^2)^2)))
		sigmax = c(sigmax,sqrt(Q[t]/4*sigmax[t]^2))
		cov = c(cov,1/2*sigmax[t]^2*(sigma^2/(sigma^2+sigmax[t]^2)+sigmax[t]^2*sigmay[t]^2/(sigma^2+sigmax[t]^2)^2+cov[t]/(sigma^2+sigmax[t]^2)))
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
		final_x_v2[i,j]=sigmax[steps]
		final_y_v2[i,j]=sigmay[steps]
		final_rho_v2[i,j] = rho[steps]
	}
}

layout(matrix(c(1:6),ncol=2,byrow=TRUE));
# xbreaks=c(0,seq(0.01,100,length.out=20))
image(sigmaxvals,sigmayvals,final_x_v2,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)));
image.scale(final_x_v2,horiz=FALSE,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))

# ybreaks=c(0,seq(0.1,10,length.out=20))
image(sigmaxvals,sigmayvals,final_y_v2,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)));
image.scale(final_y_v2,horiz=FALSE,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)))

# rhobreaks=c(0,seq(0.5,1,length.out=20))
image(sigmaxvals,sigmayvals,final_rho_v2,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)));
image.scale(final_rho_v2,horiz=FALSE,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)))

# #### until here

# ### song genetic, preference from father: sigmax2, sigmay2 image using sigmax2 not sigmax and new recursion_all function

# sigma2 = 0.25
# rho_init=0.7
# steps = 500

# sigmax2vals = seq(0,100,length.out=20)
# xx = length(sigmax2vals)
# sigmay2vals = seq(0,100,length.out=20)
# xy = length(sigmay2vals)

# final_x = array(NA,dim=c(xx,xy))
# final_y = array(NA,dim=c(xx,xy))
# final_rho = array(NA,dim=c(xx,xy))

# for(i in 1:xx){
	# for(j in 1:xy){
		# sigmax2_init = sigmax2vals[i]
		# sigmay2_init = sigmay2vals[j]		
		# variables_mat = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
		# final_x[i,j]=sqrt(variables_mat[8,1,steps])
		# final_y[i,j]=sqrt(variables_mat[8,2,steps])
		# final_rho[i,j]=variables_mat[8,4,steps]
	# }
# }

# layout(matrix(c(1:6),ncol=2,byrow=TRUE));
# xbreaks=c(0,10^seq(-30,0,length.out=20))
# image(sigmax2vals,sigmay2vals,final_x,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)));
# image.scale(final_x,horiz=FALSE,breaks=xbreaks,col=c('black',heat.colors(length(xbreaks)-2)))

# xbreaks=c(0,10^seq(-30,0,length.out=20))
# image(sigmax2vals,sigmay2vals,final_y,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)));
# image.scale(final_y,horiz=FALSE,breaks=ybreaks,col=c('black',heat.colors(length(ybreaks)-2)))

# rhobreaks=c(0,seq(0.5,0.99,length.out=20),1)
# image(sigmax2vals,sigmay2vals,final_rho,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)));
# image.scale(final_rho,horiz=FALSE,breaks=rhobreaks,col=c('black',heat.colors(length(rhobreaks)-2)))



# ### # # # song genetic, preference from mother
# ### plot contour of Q, plot slow manifold of rho / cov and sigamx2
# steps = 100
# sigma2 = 1.2
# # sigmay2 = ((30+sqrt(864))/18)*sigma2
# sigmay2 = 5.5
# poss=polyroot(c(2*sigma2^2,5*sigma2-3*sigmay2,3))

# # par(mfrow=c(2,1))

# vec=seq(0,5,length.out=100)
# covvec = seq(0,5,length.out = 100)
# Qmat = array(NA,dim=c(length(vec),length(rhovec)))
# Dmat = array(NA,dim=c(length(vec),length(rhovec)))
# testmat = array(NA,dim=c(length(vec),length(rhovec)))
# sigmay2 = sigmay2_init
# for(i in 1:length(vec)){
	# for(j in 1:length(covvec)){
		# sigmax2 = vec[i]
		# cov = covvec[j]
		# rho = cov/sqrt(sigmax2*sigmay2)
		# Q = (sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*rho*sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+1
		# Qmat[i,j] = Q
		# Dmat[i,j] = sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+(1-sqrt(Q))*rho
		# # testmat[i,j] = (sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+rho)-(Q)*sqrt(sigmax2*sigmay2)/(sigma2+Q*sigmax2)
		# testmat[i,j] =sigmax2*sigmay2/(sigma2+Q*sigmax2)
		# # testmat[i,j] = 3*sigmax2^2+sigma2*(4*sigma2-3*sigmay2+sigma2)+2*sigma2^2
		# # testmat[i,j] = (1-rho^2)*(Q)+sigma2/(sigma2+sigmax2)
		# # testmat[i,j] = (rho+sqrt(sigmax2*sigmay2)/(sigma2+sigmax2))/sqrt(Q)
	# }
# }
# contour(vec,covvec,Qmat)

# keep = c()

# for(rho in seq(0,1,length.out=10)){
# for(sigmax2_init in seq(0,10,length.out=20)){

# sigmax2=sigmax2_init
# cov = rho*sqrt(sigmay2*sigmax2)

# for(t in 1:steps){
	# sigmax2[t+1]=sigmax2[t]/4*((sigma2*(sigma2+sigmax2[t])+sigmax2[t]*sigmay2)/(sigma2+sigmax2[t])^2+2*rho[t]*sqrt(sigmax2[t]*sigmay2)/(sigma2+sigmax2[t])+1)
	# cov[t+1]=1/2*(sigmax2[t]*sigmay2/(sigma2+sigmax2[t])+rho[t]*sqrt(sigmax2[t]*sigmay2))
	# rho[t+1]=cov[t+1]/sqrt(sigmax2[t+1]*sigmay2)
# }

# Q = (sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*rho*sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+1
# start = 10
# end = 50
# lines(sigmax2[start:end],cov[start:end],col='red')
# points(sigmax2[steps],cov[steps],col='red')
# keep = rbind(keep,cbind(sigmax2[start:end],rho[start:end],cov[start:end]))
# }
# }

# # lines(vec,sqrt(vec*sigmay2)/(sigma2+vec),col='blue')
# lines(vec,(vec*sigmay2)/(sigma2+vec),col='blue')

# keep = data.frame(keep[-which(is.na(keep[,1])),])
# names(keep)=c('sigmax2','rho','cov')
# keep$inv=1/(sigma2+keep$sigmax2)
# keep$sigmax4=keep$sigmax2^2
# keep$Q = with(keep,(sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*rho*sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+1)
# keep$Qprod = keep$Q*keep$sigmax2
# keep$test = keep$sigmax2*sigmay2/(sigma2+keep$sigmax2)
# keep_ord = keep[order(keep$sigmax2),]

# # mod=lm(keep$cov/keep$Q~keep$sigmax2);
# # s=summary(mod);
# # print(s$coefficients)

# subset=which(keep$sigmax2>Re(poss[1]));mod=lm(keep$cov[subset]~keep$sigmax2[subset]);
# mod=lm(keep$cov-keep$sigmax2*sigmay2/(sigma2+keep$sigmax2)~keep$sigmax2)
# print(summary(mod)$coefficients)

# subset=which(keep$sigmax2>Re(poss[1]))
# mod=lm((keep$sigmax2[subset]+2*keep$cov[subset])~keep$sigmax2[subset])
# s=summary(mod)
# print(s$coefficients)

# mod=lm(keep$cov~keep$sigmax2+keep$Q+keep$sigmax2*keep$Q+keep$inv)
# coeffs=summary(mod)$coefficients[,1]
# print(cbind(names(summary(mod)$coefficients[,1]),(matrix(coeffs,ncol=1)))[c(1,3,4),])
# # plot(keep$sigmax2,keep$cov)
# # points(keep$sigmax2,fitted(mod),col='red')
# print(c(coeffs[1]+coeffs[4],coeffs[3]))

#####  did the last lm model above for various sigma2, to find coefficients as a function of sigma2, which only resulted in finding a way to rewrite cov in terms of sigmax2 and Q
# steps = 100
# # sigma2 = 1.1
# # sigmay2_init = ((30+sqrt(864))/18)*sigma2
# sigmay2_init = 5
# sigma2_vals = seq(1,1.3,length.out=50)
# coeffs=array(NA,c(5,length(sigma2_vals)))
# for(i in 1:length(sigma2_vals)){
# sigma2 = sigma2_vals[i]
# keep = c()

# for(rho in seq(0,1,length.out=10)){
# for(sigmax2_init in seq(0,10,length.out=20)){
	# sigmay2 = sigmay2_init
	
	# sigmax2=sigmax2_init
	# cov = rho*sqrt(sigmay2*sigmax2)
	
	# for(t in 1:steps){
		# sigmax2[t+1]=sigmax2[t]/4*((sigma2*(sigma2+sigmax2[t])+sigmax2[t]*sigmay2)/(sigma2+sigmax2[t])^2+2*rho[t]*sqrt(sigmax2[t]*sigmay2)/(sigma2+sigmax2[t])+1)
		# cov[t+1]=1/2*(sigmax2[t]*sigmay2/(sigma2+sigmax2[t])+rho[t]*sqrt(sigmax2[t]*sigmay2))
		# rho[t+1]=cov[t+1]/sqrt(sigmax2[t+1]*sigmay2)
	# }
	
	# Q = (sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*rho*sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+1
	
	# start = 10
	# end = 50
	# keep = rbind(keep,cbind(sigmax2[start:end],rho[start:end],cov[start:end]))
# }
# }

# keep = data.frame(keep[-which(is.na(keep[,1])),])
# names(keep)=c('sigmax2','rho','cov')
# keep$inv=1/(sigma2+keep$sigmax2)
# keep$Q = with(keep,(sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*rho*sqrt(sigmax2*sigmay2)/(sigma2+sigmax2)+1)
# keep_ord = keep[order(keep$sigmax2),]

# mod=lm(keep$cov~keep$sigmax2+keep$Q+keep$sigmax2*keep$Q+keep$inv)
# coeffs[,i]=summary(mod)$coefficients[,1]
# }

