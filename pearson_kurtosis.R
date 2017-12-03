setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
source('saveit.R')
source('init_conds.R')
library(RColorBrewer)
library(pracma)
library(PearsonDS)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 193
sigmay2 = 2
sigmax2 = 0.8
sigma2 = 1.1
mut_prob = 0.0
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
minprob = 0

steps = 30000
store = steps+1
source('range_setup_long.R')

excess_kurt_vals = c(-0.2,-0.15,-0.1,-0.05)
mut_prob_vals = c(0,0.001)
xk = length(excess_kurt_vals)
xm = length(mut_prob_vals)


ic=init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,NA,NA)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

# p1 = dynamics_memory()

# equilibrium = as.list(1:(xk*xm))
# dim(equilibrium)=c(xk,xm)

# for(k in 1:xk){
	# for(m in 1:xm){
	# excess_kurt = excess_kurt_vals[k]
	# mut_prob = mut_prob_vals[m]
	
	# f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
	# f_init2 = f_init/sum(f_init)
	
	# f_init = f_init2
	# p2 = dynamics_memory()
	
	# equilibrium[[k,m]] = p2$Pm
	# }
# }


# saveit(p1=p1,equilibrium=equilibrium,excess_kurt_vals=excess_kurt_vals,mut_prob_vals=mut_prob_vals,sigmax2=sigmax2,sigma2=sigma2,sigmay2=sigmay2,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/pearson_kurtosis.Rdata')

#### 
load('/Users/eleanorbrush/Documents/research/song_learning_evolution/pearson_kurtosis.Rdata')

col_vec = c('black',brewer.pal(9,'Set1')[1])
lwd = 2
marg = c(0.53,0.4,0.0,0.15)
omarg = c(0.03,1,0.35,0.0)
ylim = c(0,0.22)

width = 6.5
height = 4.5

w = 97+(-27:27)

pdf('/Users/eleanorbrush/Desktop/pearson_kurtosis.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,1,0))
layout(matrix(1:4,ncol=2,byrow=TRUE))

for(i in 1:xk){
	plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
	for(j in 1:2){
		lines(mrange[w]+1,equilibrium[[i,j]][w],lwd=lwd,col=col_vec[j])
	}
	mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
	mtext('Frequency',side=2,line=1.9,at=mean(ylim),cex=largefontsize/smallfontsize)
	if(i ==4){		legend(-7.5,0.23,legend=mut_prob_vals[1:2],lty=1,col=col_vec[1:2],lwd=lwd,bty='n')}
}

dev.off()