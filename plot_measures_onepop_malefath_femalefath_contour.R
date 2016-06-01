library(ggplot2)
library(RColorBrewer)
library(reshape)
library(gridExtra)
source('local_max.R')
source('entropy.R')
source('ind2sub.R')
source('int.R')
source('get_legend.R')

mypal=brewer.pal(8,'Set1')

d = dim(Pm_onepop)
P = prod(d)
Tend=d[2]


Ns = length(sigma2_vals)
Nfs = length(fmix_sigma2_vals)
Nms = length(mmix_sigma2_vals)
Nmp = length(mut_prob_vals)

subset = 601:1201;
subset = 1:Nm;
Tend = dim(Pm_onepop[[1]])[2]

thresh = 1e-2
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
	ex = int(mrange*Pm_onepop[[i]][,Tend])
	vx = int((mrange-(ex))^2*Pm_onepop[[i]][,Tend])
	var_mat_m[i] = vx
	ex = int(mrange*Pf_onepop[[i]][,Tend])
	vx = int((mrange-(ex))^2*Pf_onepop[[i]][,Tend])
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
	smat[i] = sigma2_vals[sub[1]]
	fmat[i] = fmix_sigma2_vals[sub[2]]
	mmat[i] = mmix_sigma2_vals[sub[3]]
	pmat[i] = mut_prob_vals[sub[4]]
	colmat[i] = mypal[sub[2]]
}

p=2
mvals=c(2,3)
lm = length(mvals)
fvals = 1:Nfs

subset = (m0-250):(m0+250)
ylim = c(0,3)
marg = c(0.25,0.1,0.25,0.1)

examples = list()
contour_sd = list()
contour_num = list()

# ms_toplot = cbind(rep(mvals,each=2),c(1,4,4,Ns))
s_toplot = c(1,3,Ns)
for(i in 1:length(s_toplot)){
	s = s_toplot[i]
	dist = c()
	m_init = dnorm(mrange[subset],mmin,mmix_sigma2_vals[m])
	m_init[m_init>ylim[2]] = ylim[2]
	dist = c(dist,m_init)
	for(m in mvals){
		for(f in fvals){
			toadd =Pm_onepop[[s,f,m,p]][subset,Tend]
			toadd[toadd>ylim[2]] = ylim[2]
			dist = c(dist,toadd)
		}
	}
	dist = data.frame(dist = dist, mrange = rep(mrange[subset]+1,length(fvals)*length(mvals)+1),f_sigma2 = as.factor(rep(c(fmix_sigma2_vals[Nfs]-1,rep(fmix_sigma2_vals[fvals],times=length(mvals))),each=length(subset))),m_sigma2=as.factor(rep(c(0,rep(mmix_sigma2_vals[mvals],each=length(fvals))),each=length(subset))))
	
examples[[i]] <- ggplot(dist,aes(x=mrange,y=dist,color=f_sigma2, linetype=m_sigma2)) + geom_line() + 
	theme_bw() +
	theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position=c(mrange[subset[length(subset)]]-1,0.58)) + 
	scale_color_manual(values=c('white',mypal[fvals]),breaks=levels(dist$f_sigma2)[2:(Nfs+1)])+
	scale_linetype_manual(values=c(1,1,2),breaks=mmix_sigma2_vals[mvals]) +
	scale_y_continuous(limits=ylim)  + xlab('Song') + ylab('Probability') + labs(color=expression(sigma[f]),linetype=expression(sigma[m])) 
	}


legend_ex <- get_legend(examples[[1]])
examples[[1]] = examples[[1]] + theme(legend.position='none')
examples[[2]] = examples[[2]] + theme(legend.position='none')
examples[[3]] = examples[[3]] + theme(legend.position='none')


xbreaks = seq(0,1,by=0.2)
ybreaks = seq(0,0.1,by=0.02)


blues=brewer.pal(9, "Blues")[3:9]
bluePalette = colorRampPalette(blues)
sdbreaks=c(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.5,1)
sdpal = bluePalette(length(sdbreaks))
reds=brewer.pal(9, "Reds")[3:9]
redPalette = colorRampPalette(reds)
numbreaks = c(2,15,30,40)
# numbreaks = 1:max(num_peaks)
numpal = redPalette(length(numbreaks))

# for(i in 1:lm){
	i=1
	m = mvals[i]
	
	toplot = melt((sqrt(var_mat_m[,1:Nfs,m,p])),varnames = c('sigma2','f_sigma2'))
	toplot$sigma2 = (rep(sigma2_vals,times=(Nfs)))
	toplot$f_sigma2 =(rep(fmix_sigma2_vals[1:Nfs],each=Ns))
	colnames(toplot)[3] = 'sd'
	toplot$num = as.vector(num_peaks[,1:Nfs,m,p])
	
	
contour_sd <- ggplot(toplot, aes(x=f_sigma2,y= sigma2, z=sd) )+ 
	stat_contour(aes(x=f_sigma2,y=sigma2,z=sd, colour=as.factor(..level..)),breaks=sdbreaks) + 
	scale_colour_manual(breaks=sdbreaks,values=sdpal) +
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position=c(0.6,.6)) + 
		 scale_x_continuous(limits=range(toplot$f_sigma2),breaks=xbreaks,labels=xbreaks)+
		scale_y_continuous(limits=c(0,0.1),breaks=ybreaks,labels=ybreaks) +xlab('Female stand. dev.') + ylab('Promiscuity')+
		labs(color='Stand. dev.')
		
contour_num <- ggplot(toplot, aes(x=f_sigma2,y= sigma2, z=num) )+ 
	stat_contour(aes(x=f_sigma2,y=sigma2,z=num, colour=as.factor(..level..)),breaks=numbreaks) + 
	scale_colour_manual(breaks=numbreaks,values=numpal) +
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position=c(0.6,.8)) + 
		 scale_x_continuous(limits=range(toplot$f_sigma2),breaks=xbreaks,labels=xbreaks)+
		scale_y_continuous(limits=c(0,0.1),breaks=ybreaks,labels=ybreaks) + xlab('Female stand. dev.') + ylab('Promiscuity') + 
		labs(color='# of peaks')

# }

legend_sd <- get_legend(contour_sd)
contour_sd = contour_sd + theme(legend.position='none') 
legend_num <- get_legend(contour_num)
contour_num = contour_num + theme(legend.position='none') 

examples_tot = arrangeGrob(examples[[1]],examples[[2]],examples[[3]],legend_ex,nrow=1,widths=c(1,1,1,0.3))
contour_tot = arrangeGrob(contour_sd,legend_sd,contour_num,legend_num,nrow=1,widths=c(1,0.3,1,0.3))

	
pdf(file=paste('/Users/eleanorbrush/Documents/research/song_learning_evolution/examples_and_summary_malefath_femalefath_cont_p=',p,'.pdf',sep=''),width=6.83,height=5)
grid.arrange(examples_tot,contour_tot,ncol=1)
dev.off()
