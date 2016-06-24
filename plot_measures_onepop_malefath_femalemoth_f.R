library(ggplot2)
library(RColorBrewer)
library(reshape)
library(gridExtra)
source('local_max.R')
source('entropy.R')
source('ind2sub.R')
source('int.R')
source('get_legend.R')

eightpal=brewer.pal(8,'Set1')

d = dim(Pm_onepop)
P = prod(d)
Tend=d[2]


Ns = length(sigma_vals)
Nfs = length(f_sigma_vals)
Nms = length(m_sigma_vals)
Nmp = length(mut_prob_vals)

subset = 601:1201;
subset = 1:Nm;
Tend = dim(Pm_onepop[[1]])[2]

thresh = 5e-4
half = floor(Nm/2)

four_freq = array(NA,dim(Pm_onepop))
four_mat = list()
num_peaks = array(NA,dim(Pm_onepop))
var_mat_m = array(NA,c(dim(Pm_onepop)))
var_mat_f = array(NA,c(dim(Pf_onepop)))
ent_mat = array(Inf,dim=dim(Pm_onepop))



for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tend]
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourCoeffs[1] / length(subset)	
	mf = mean(coeffs[1:200])
	ml = intersect(local_max(round(coeffs[1:half],5)),which(coeffs[1:half]>mf/2))
	if(length(ml)==1 || length(ml)==0){four_freq[i] = 1} else{ four_freq[i] = median(diff(ml))}
	four_mat[[i]] = coeffs[1:half]
	m = local_max(eq_pop)
	m = intersect(m, which(eq_pop>thresh))
	l = length(m)
	num_peaks[i] = l
	sub = ind2sub(dim(Pm_onepop),i)
	ex = sum(mrange*Pm_onepop[[i]][,Tend])
	vx = sum((mrange-(ex))^2*Pm_onepop[[i]][,Tend])
	var_mat_m[i] = vx
	ex = sum(mrange*Pf_onepop[[i]][,Tend])
	vx = sum((mrange-(ex))^2*Pf_onepop[[i]][,Tend])
	var_mat_f[i] = vx
	ent_mat[i] = ent(eq_pop)
}
dim(four_mat) = dim(Pm_onepop)

smat= array(0,dim=dim(Pm_onepop))
fmat= array(0,dim=dim(Pm_onepop))
mmat= array(0,dim=dim(Pm_onepop))
pmat = array(0,dim=dim(Pm_onepop))
colmat = array(0,dim=dim(Pm_onepop))

for(i in 1:P){
	sub = ind2sub(dim(Pm_onepop),i)
	smat[i] = sigma_vals[sub[1]]
	fmat[i] = f_sigma_vals[sub[2]]
	mmat[i] = m_sigma_vals[sub[3]]
	pmat[i] = mut_prob_vals[sub[4]]
	colmat[i] = eightpal[sub[2]]
}

p=2
mvals=c(1,2)
lm = length(mvals)
fvals = c(1:3)
svals = c(1,3,5,6)

subset = (m0-300):(m0+300)
ylim = c(0,0.05)
marg = c(0.25,0.1,0.25,0.1)

examples = list()
bubble = list()

mf_toplot = cbind(rep(mvals,each=2),c(1,Nfs,2,Nfs))
for(i in 1:dim(mf_toplot)[1]){
	m = mf_toplot[i,1]
	f = mf_toplot[i,2]
	dist = c()
	m_init = dnorm(mrange[subset],mmin,m_sigma_vals[m])
	m_init[m_init>ylim[2]] = ylim[2]
	dist = c(dist,m_init)
	for(s in svals){
		toadd =Pm_onepop[[s,f,m,p]][subset,Tend]
		toadd[toadd>ylim[2]] = ylim[2]
		dist = c(dist,toadd)
	}
	dist = data.frame(dist = dist, mrange = rep(mrange[subset]+1,length(svals)+1),sigma = as.factor(rep(c(0,sigma_vals[svals]),each=length(subset))))
	
examples[[i]] <- ggplot(dist,aes(x=mrange,y=dist,color=sigma)) + geom_line() + 
	theme_bw() +
	theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position='none') + 
	scale_color_manual(values=c('black',eightpal[1:length(svals)]))+
	scale_y_continuous(limits=ylim)  + xlab('') + ylab('')
}

examples[[3]] = examples[[3]] + xlab('Song') + ylab('Probability')
examples[[4]] = examples[[4]] + xlab('Song') + ylab('Probability')

eightpal2= eightpal
eightpal2[svals] = eightpal[1:length(svals)]
eightpal2[setdiff(1:length(eightpal),svals)] = eightpal[(length(svals)+1):length(eightpal)]

for(i in 1:lm){
	m = mvals[i]
	
	toplot = melt(log(sqrt(var_mat_m[(1:Ns),1:Nfs,m,p]),base=10),varnames = c('sigma','f_sigma'))
	toplot$sigma = as.factor(rep((sigma_vals),times=(Nfs)))
	toplot$f_sigma = log(rep(f_sigma_vals[1:Nfs],each=Ns),base=10)
	colnames(toplot)[3] = 'y'
	toplot$size = as.vector(num_peaks[(1:Ns),1:Nfs,m,p])

	xbreaks = 10^(-3:0)
	ybreaks = 10^(-10:0)
	
	bubble[[i]] <- ggplot(toplot, aes(x=f_sigma, y=y, size=size, color = sigma),guide=FALSE)+
		# geom_hline(yintercept = log(f_sigma_vals[1:Nfs],base=10), color =eightpal[1:(Nfs)],size=.2) + 
		geom_line(size=0.5) + geom_point()+ scale_size_area(max_size = 10,limits = c(0,max(num_peaks)))+
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm")) + 
		scale_color_manual(values=(eightpal2))+
		labs(color=expression(sigma),size='# of peaks') +
		 scale_x_continuous(limits=range(log(f_sigma_vals[1:Nfs],base=10))+c(-.4,.1),breaks=log(xbreaks,base=10),labels=xbreaks)+
		scale_y_continuous(limits=c(floor(min(toplot$y)),0.5),breaks=log(ybreaks,base=10),labels=ybreaks)+
		 guides(size = FALSE) 
}


bubble[[1]] = bubble[[1]] + theme(legend.position='none') + xlab('') + ylab('')
legend_bub <- get_legend(bubble[[2]])
bubble[[2]] = bubble[[2]] + theme(legend.position='none') + xlab('Female standard deviation') + ylab('Male standard deviation')

	
# pdf(file=paste('/Users/eleanorbrush/Documents/research/song_learning_evolution/examples_and_summary_malefath_femalemoth_f_p=',p,'.pdf',sep=''),width=6.83,height=5)
grid.arrange(examples[[1]],examples[[2]],bubble[[1]],legend_bub,examples[[3]],examples[[4]],bubble[[2]],ncol=4,widths=c(1,1,1,.4))
# dev.off()

