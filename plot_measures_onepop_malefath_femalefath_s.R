# source('plot_setup.R')

Ncol = 8 
mypal=(brewer.pal(Ncol,'Set1'))
mypal=brewer.pal(Ncol,'Spectral')


p=1

mvals=c(1,2)
lm = length(mvals)
fvals = apply(matrix(c(0.001,1),nrow=1),2,function(x) which(is.element(f_sigma_vals,x)))

subset = (m0-250):(m0+250)
ylim = c(0,0.05)
marg = c(0.25,0.1,0.25,0.1)

examples = list()
bubble = list()

ms_toplot = cbind(rep(mvals,each=3),rep(c(1,2,Ns),2))
for(i in 1:dim(ms_toplot)[1]){
	m = ms_toplot[i,1]
	s = ms_toplot[i,2]
	dist = c()
	m_init = dnorm(mrange[subset],mmin,m_sigma_vals[m])
	m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
	m_init = m_init / sum(m_init)
	m_init[m_init>ylim[2]] = ylim[2]
	dist = c(dist,m_init)
	for(f in rev(fvals)){
		toadd =Pm_onepop[[s,f,m,p]][subset,Tend]
		toadd[toadd>ylim[2]] = ylim[2]
		dist = c(dist,toadd)
	}
	dist = data.frame(dist = dist, mrange = rep(mrange[subset]+1,length(fvals)+1),f_sigma = as.factor(rep(c(f_sigma_vals[Nfs]-1,f_sigma_vals[fvals]),each=length(subset))))
	
examples[[i]] <- ggplot(dist,aes(x=mrange,y=dist,color=f_sigma)) + geom_line() + 
	theme_bw() +
	theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position='none') + 
	scale_color_manual(values=c('white',rev(mypal[fvals])))+
	scale_y_continuous(limits=ylim)  + xlab('') + ylab('')
}

examples[[4]] = examples[[4]] + xlab('Song') + ylab('Probability')
examples[[5]] = examples[[5]] + xlab('Song') + ylab('Probability')
examples[[6]] = examples[[6]] + xlab('Song') + ylab('Probability')

for(i in 1:lm){
	m = mvals[i]
	
	toplot = melt(log(sqrt(var_mat_m[,1:Nfs,m,p]),base=10),varnames = c('sigma','f_sigma'))
	toplot$sigma = log(rep(sigma_vals,times=(Nfs)),base=10)
	toplot$f_sigma = as.factor(rep(f_sigma_vals[1:Nfs],each=Ns))
	colnames(toplot)[3] = 'y'
	toplot$size = as.vector(num_peaks[,1:Nfs,m,p])
	
	xbreaks = 10^(-3:0)
	ybreaks = 10^(-10:0)
	
	bubble[[i]] <- ggplot(toplot, aes(x=sigma, y=y, size=size, color = f_sigma),guide=FALSE)+
		geom_hline(yintercept = log(f_sigma_vals[1:Nfs],base=10), color =mypal[1:(Nfs)],size=.2) + 
		geom_line(size=0.5) + geom_point()+ scale_size_area(max_size = 5,limits = c(0,max(num_peaks)))+
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm")) + 
		scale_color_manual(values=mypal)+
		labs(color=expression(sigma[f]),size='# of peaks') +
		 scale_x_continuous(limits=range(log(sigma_vals,base=10))+c(-.4,.1),breaks=log(xbreaks,base=10),labels=xbreaks)+
		scale_y_continuous(limits=c(floor(min(toplot$y)),0.5),breaks=log(ybreaks,base=10),labels=ybreaks)+
		 guides(size = FALSE) 
}


bubble[[1]] = bubble[[1]] + theme(legend.position='none') + xlab('') + ylab('')
legend_bub <- get_legend(bubble[[2]])
bubble[[2]] = bubble[[2]] + theme(legend.position='none') + xlab('Promiscuity') + ylab('Standard deviation')

	
pdf(file=paste('/Users/eleanorbrush/Documents/research/song_learning_evolution/examples_and_summary_malefath_femalefath_s_p=',p,'.pdf',sep=''),width=6.83,height=5)
grid.arrange(examples[[1]],examples[[2]],examples[[3]],bubble[[1]],legend_bub,examples[[4]],examples[[5]],examples[[6]],bubble[[2]],ncol=5,widths=c(1,1,1,1,.4))
dev.off()

