setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_by_mode.R')
source('saveit.R')
library(RColorBrewer)

trait_chunk_num = 75
source('range_setup.R')
mrange=c(rev(seq(0,-10,by=-0.25)),c(seq(0.25,8,by=0.25)))
Nm = length(mrange)
midpt = which(mrange==-1)
frange=mrange
Nf = Nm
sigmay2 = 2
sigmax2 = 0.1
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob = 0.01

sigma2 = 1.5

k1 = 9
k2 = 12

ic = init_conds('norm','norm','step',sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight1 = ic$fixed_weight
p = unique(fixed_weight1)

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,NA,NA)
fixed_weight2 = ic$fixed_weight


w = setdiff(which(round(mrange,3)>=-7),which(round(mrange,3)>5))

fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12
col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.6,0.5,0.02,0.15)
omarg = c(0.1,0.4,0.6,0.0)

width = 6.8
height = 3.25

pdf(file='/Users/eleanorbrush/Desktop/step_function_example.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.5,0))

plot(mrange[w]+1,fixed_weight2[w],t='o',xlim=range(mrange[w]+1)+c(-.01,.01),yaxt='n',xlab='',ylab='',lwd=lwd,col=col_vec[1],ylim=range(c(fixed_weight1,fixed_weight2))+c(-.00,0),xaxt='n')
points(mrange[w]+1,fixed_weight1[w],t='o',col=col_vec[2],lwd=lwd)
axis(2,at=c(p[2],p[3],c(0,0.04,0.08)),labels=c(expression(p[1]),expression(p[2]),c(0,0.04,0.08)))
axis(1,at=-6:6,labels=-6:6)
mtext('Difference in preference and song, y-x',side=1,line=1.6,at=-0.3,cex=largefontsize/smallfontsize)
mtext('Preference',side=2,line=1.7,at=p[3]/2,cex=largefontsize/smallfontsize)

dev.off()