# setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
source('saveit.R')
source('recursion_all.R')
source('init_conds.R')
library(RColorBrewer)
library(pracma)
library(PearsonDS)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 281
minweight = 10^(-320)
mut_prob_nonzero = 0.001
source('range_setup.R')

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
rho = 0

steps = 30000
store_vec = c(1,100,10000,20000,30000)

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,NA,NA)
m_init_norm = ic$m_init
f_init_norm = ic$f_init
fixed_weight = ic$fixed_weight	

k1 = 21
k2 = 15

ic = init_conds('step','norm','norm',sigmax2,sigmay2,sigma2,k1,k2)
m_init_step = ic$m_init


k1 = 21
k2 = 15

ic = init_conds('norm','step','norm',sigmax2,sigmay2,sigma2,k1,k2)
f_init_step = ic$f_init

f_init = f_init_norm
m_init = m_init_norm
mut_prob = mut_prob_nonzero
p2 = dynamics_memory(store=TRUE)

f_init = f_init_norm
m_init = m_init_step

mut_prob = 0
p1_step_song_dist = dynamics_memory(store=TRUE)	

mut_prob = mut_prob_nonzero
p2_step_song_dist = dynamics_memory(store=TRUE)

f_init = f_init_step
m_init = m_init_norm

mut_prob = 0
p1_step_pref_dist = dynamics_memory(store=TRUE)	

mut_prob = mut_prob_nonzero
p2_step_pref_dist = dynamics_memory(store=TRUE)

saveit(p2=p2,p1_step_song_dist=p1_step_song_dist,p2_step_song_dist=p2_step_song_dist,p1_step_pref_dist=p1_step_pref_dist,p2_step_pref_dist=p2_step_pref_dist,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,store_vec=store_vec,trait_chunk_num=trait_chunk_num,file='/homes/ebrush/priv/song_learning_evolution/mutation_test.Rdata')

# load('/Users/eleanorbrush/Documents/research/song_learning_evolution/mutation_test.Rdata')

# r = recursion_all(sigmax2,sigmay2,sigma2,rho)

# t_toplot = c(1:4)

# p = c()

# for(t in store_vec[t_toplot]){
	# p = cbind(p,dnorm(mrange,mean=-1,sd=sqrt(r[3,1,t]))/sum(dnorm(mrange,mean=-1,sd=sqrt(r[3,1,t]))))
# }

# col_vec = brewer.pal(11,'RdBu')[c(11,9,3,2)]
# lwd = 2
# marg = c(0.42,0.3,0.0,0.15)
# omarg = c(0.03,1,1.7,0.0)

# width = 6.5
# height = 5
# ylim = c(0,0.2)
# xlim=c(-6,6)
# t = dim(p2$Pm_store)[2]
# w = which(p2$Pm_store[,t]>1e-10)

# pdf('/Users/eleanorbrush/Desktop/mutation_sensitivity.pdf',width=width,height=height,family=fontfamily)

# par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
# layout(matrix(1:6,ncol=3,byrow=FALSE))

# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p[w,t],lwd=lwd,col=col_vec[t])
# }
# # mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)
# mtext('Both traits normal',side=3,line=0.1,cex=largefontsize/smallfontsize)
# legend(-4.5,0.1,legend=store_vec[t_toplot],lty=1,col=col_vec,lwd=lwd,bty='n')

# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p2$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
# }
# mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p1_step_song_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
# }
# # mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# # mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)
# mtext('Songs follow step function',side=3,line=0.1,cex=largefontsize/smallfontsize)

# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p2_step_song_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
# }
# mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# # mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)


# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p1_step_pref_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
# }
# # mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# # mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)
# mtext("Pref.'s follow step function",side=3,line=0.1,cex=largefontsize/smallfontsize)

# plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='',xlim=xlim)
# for(t in 1:length(t_toplot)){
	# lines(mrange[w]+1,p2_step_pref_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
# }
# mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
# # mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

# dev.off()