source('dynamics_pxy.R')
source('dynamics_by_mode_new_numbers.R')
source('init_conds.R')
source('saveit.R')
library(RColorBrewer)
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 281
minweight = 10^(-320)
mut_prob = 0
source('range_setup.R')

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
mut_prob = 0.0
rho = 0

store = 1

k1 = 21
k2 = 15

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,k1,k2)
f_init = ic$f_init
fixed_weight = ic$fixed_weight
m_init1 = ic$m_init

ic = init_conds('step','norm','norm',sigmax2,sigmay2,sigma2,k1,k2)
m_init2 = ic$m_init

# m_init = m_init1
# steps = 1
# p1_hold = dynamics_pxy()
# steps = 30000
# p1 = dynamics_mode3()
# p1$pxy = p1_hold$pxy
# p1$z = p1_hold$z

# m_init = m_init2
# steps = 1
# p2_hold = dynamics_pxy()
# steps = 30000
# p2 = dynamics_mode3()
# p2$pxy = p2_hold$pxy
# p2$z = p2_hold$z

# saveit(p1=p1,p2=p2,k1=k1,k2=k2,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,file='/Users/eleanorbrush/Desktop/peak_example_step_song_dist.Rdata')

# load('/Users/eleanorbrush/Documents/research/song_learning_evolution/peak_example_step_song_dist.Rdata')

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.45,0.43,0.02,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 4.5

pdf('/Users/eleanorbrush/Desktop/peak_example_step_song_dist.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:4,ncol=2,byrow=TRUE))

t=1
w1 = which(p1$Pm[,t]>1e-10)
w2 = which(p2$Pm[,t]>1e-10)
w1 = which(abs(mrange+1)<6.1)
w2 = which(abs(mrange+1)<6.1)

plot(mrange[w1]+1,m_init1[w1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,0.04),yaxt='n')
axis(2,at=seq(0,0.1,0.02))
lines(mrange[w1]+1,m_init2[w1],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Song, x',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste("Frequency, ",P[m] ,'(x)')),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p1$Pm[w1,steps],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,0.35),yaxt='n')
axis(2,at=seq(0,0.4,by=0.05))
points(mrange+1,p2$Pm[,steps],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Song, x',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste("Frequency, ",P[m] ,'(x)')),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p1$z[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,0.033))
points(mrange[w1]+1,p2$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Preference, y',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext((expression(paste('Frac. of attractive males, ',Z[y]))),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,log(p1$Pf[w1,t]/p1$z[w1,t]),t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',yaxt='n',ylim=log(c(0.7,4.4)))
points(mrange[w1]+1,log(p2$Pf[w1,t]/p2$z[w1,t]),lwd=lwd,col='black',t='l')
axis(2,at=log(2^(-3:3)),labels=2^(-3:3))
mtext('Preference, y',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste(P[f] ,'(y)',' / ',Z[y])),side=2,line=1.7,cex=largefontsize/smallfontsize)

dev.off()

####

k1 = 21
k2 = 15

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
fixed_weight = ic$fixed_weight
f_init1 = ic$f_init

ic = init_conds('norm','step','norm',sigmax2,sigmay2,sigma2,k1,k2)
f_init2 = ic$f_init

# f_init = f_init1
# steps = 1
# p1_hold = dynamics_pxy()
# steps = 30000
# p1 = dynamics_mode3()
# p1$pxy = p1_hold$pxy
# p1$z = p1_hold$z

# f_init = f_init2
# steps = 1
# p2_hold = dynamics_pxy()
# steps = 30000
# p2 = dynamics_mode3()
# p2$pxy = p2_hold$pxy
# p2$z = p2_hold$z

# saveit(p1=p1,p2=p2,k1=k1,k2=k2,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,file='/Users/eleanorbrush/Desktop/peak_example_step_pref_dist.Rdata')

# load('/Users/eleanorbrush/Documents/research/song_learning_evolution/peak_example_step_pref_dist.Rdata')

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.45,0.43,0.02,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 4.5

pdf('/Users/eleanorbrush/Desktop/peak_example_step_pref_dist.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:4,ncol=2,byrow=TRUE))

t=1
w1 = which(p1$Pm[,t]>1e-15)
w2 = which(p2$Pm[,t]>1e-15)
w1 = which(abs(mrange+1)<6.1)
w2 = which(abs(mrange+1)<6.1)

plot(mrange[w1]+1,f_init1[w1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,0.04),yaxt='n')
axis(2,at=seq(0,0.1,0.02))
points(mrange[w1]+1,f_init2[w1],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Preference, y',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste("Frequency, ",P[f] ,'(y)')),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p1$Pm[w1,steps],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,0.35),yaxt='n')
axis(2,at=seq(0,0.4,by=0.05))
points(mrange+1,p2$Pm[,steps],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Song, x',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste("Frequency, ",P[m] ,'(x)')),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p2$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='',ylim=range(c(p1$z[w1,t],p2$z[w2,t])))
# points(mrange[w1]+1,p1$z[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='')
mtext('Preference, y',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext((expression(paste('Frac. of attractive males, ',Z[y]))),side=2,line=1.7,cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,log(p1$Pf[w1,t]/p1$z[w1,t]),t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',yaxt='n',ylim=log(c(0.7,4.4)))
lines(mrange[w1]+1,log(p2$Pf[w1,t]/p2$z[w1,t]),lwd=lwd,col='black')
axis(2,at=log(2^(-3:3)),labels=2^(-3:3))
mtext('Preference, y',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext(expression(paste(P[f] ,'(y)',' / ',Z[y])),side=2,line=1.7,cex=largefontsize/smallfontsize)

dev.off()

