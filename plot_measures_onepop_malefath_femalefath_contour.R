# source('plot_setup.R')
Ncol = 8 
mypal=(brewer.pal(Ncol,'Set1'))
mypal=brewer.pal(Ncol,'Spectral')

p=2
mvals=c(1,2)
lm = length(mvals)
fvals = apply(matrix(c(0.001,1),nrow=1),2,function(x) which(is.element(f_sigma_vals,x)))

subset = (m0-300):(m0+300)
ylim = c(0,0.01)
marg = c(0.25,0.1,0.25,0.1)

examples = list()
contour_sd = list()
contour_num = list()

# ms_toplot = cbind(rep(mvals,each=2),c(1,4,4,Ns))
s_toplot = c(1,Ns)
for(i in 1:length(s_toplot)){
	s = s_toplot[i]
	dist = c()
	m_init = dnorm(mrange[subset],mmin,m_sigma_vals[m])
	m_init[m_init>ylim[2]] = ylim[2]
	dist = c(dist,m_init)
	for(m in mvals){
		for(f in rev(fvals)){
			toadd =Pm_onepop[[s,f,m,p]][subset,Tend]
			toadd[toadd>ylim[2]] = ylim[2]
			dist = c(dist,toadd)
		}
	}
	dist = data.frame(dist = dist, mrange = rep(mrange[subset]+1,length(fvals)*length(mvals)+1),f_sigma = as.factor(rep(c(f_sigma_vals[Nfs]-1,rep(f_sigma_vals[fvals],times=length(mvals))),each=length(subset))),m_sigma=as.factor(rep(c(0,rep(m_sigma_vals[mvals],each=length(fvals))),each=length(subset))))
	
examples[[i]] <- ggplot(dist,aes(x=mrange,y=dist,color=f_sigma, linetype=m_sigma)) + geom_line() + 
	theme_bw() +
	theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"))+#,legend.position=c(mrange[subset[length(subset)]]-0.5,0.58)) + 
	scale_color_manual(values=c('white',rev(mypal[fvals])),breaks=levels(dist$f_sigma)[2:(Nfs+1)])+
	scale_linetype_manual(values=c(1,1,0),breaks=m_sigma_vals[mvals]) +
	scale_y_continuous(limits=ylim)  + xlab('Song') + ylab('Probability') + labs(color=expression(sigma[f]),linetype=expression(sigma[m]))+guides(linetype = FALSE)  }
	
# i=2
# s = s_toplot[i]
	# dist = c()
	# m_init = dnorm(mrange[subset],mmin,m_sigma_vals[m])
	# m_init[m_init>ylim[2]] = ylim[2]
	# dist = c(dist,m_init)
	# for(m in mvals){
		# for(f in fvals){
			# toadd =Pm_onepop[[s,f,m,p]][subset,Tend]
			# toadd[toadd>ylim[2]] = ylim[2]
			# dist = c(dist,toadd)
		# }
	# }
	# dist = data.frame(dist = dist, mrange = rep(mrange[subset]+1,length(fvals)*length(mvals)+1),f_sigma = as.factor(rep(c(f_sigma_vals[Nfs]-1,rep(f_sigma_vals[fvals],times=length(mvals))),each=length(subset))),m_sigma=as.factor(rep(c(0,rep(m_sigma_vals[mvals],each=length(fvals))),each=length(subset))))
	
# examples[[i]] <- ggplot(dist,aes(x=mrange,y=dist,color=f_sigma, linetype=m_sigma)) + geom_line() + 
	# theme_bw() +
	# theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"))+#,legend.position=c(mrange[subset[length(subset)]]-0.5,0.58)) + 
	# scale_color_manual(values=c('white',mypal[fvals]),breaks=levels(dist$f_sigma)[2:(Nfs+1)])+
	# scale_linetype_manual(values=c(1,1,0),breaks=m_sigma_vals[mvals]) +
	# scale_y_continuous(limits=ylim)  + xlab('Song') + ylab('Probability') + labs(color=expression(sigma[f]),linetype=expression(sigma[m]))

legend_ex <- get_legend(examples[[2]])
examples[[1]] = examples[[1]] + theme(legend.position='none')
examples[[2]] = examples[[2]] + theme(legend.position='none')
# examples[[3]] = examples[[3]] + theme(legend.position='none')

examples_tot = arrangeGrob(examples[[1]],examples[[2]],legend_ex,nrow=1,widths=c(1,1,0.3))

xbreaks = seq(0,1,by=0.2)
ybreaks = seq(0,2,by=0.1)
ylim=c(0,0.5)

blues=brewer.pal(9, "Blues")[3:9]
bluePalette = colorRampPalette(blues)
grad = seq(0,1,by=0.01)
sdbreaks=c(seq(0,1,by=0.1))
sdpal_base = bluePalette(length(grad))
sdpal = sdpal_base[apply(matrix(sdbreaks,nrow=1),2,function(x) which(round(grad,2) == round(x,2)))]
sdbreaks_toshow = seq(0,1,by=0.1)
sdpal_toshow = sdpal_base[apply(matrix(sdbreaks_toshow,nrow=1),2,function(x) which(round(grad,2) == round(x,2)))]
reds=brewer.pal(9, "Reds")[3:9]
redPalette = colorRampPalette(reds)
numbreaks = c(2,10,seq(20,100,by=20))
# numbreaks = 1:max(num_peaks)
set = min(numbreaks):max(numbreaks)
numpal = redPalette(length(set))
numpal = numpal[apply(matrix(numbreaks,nrow=1),2, function(x) which(set == round(x,2)))]

# for(i in 1:lm){
	i=1
	m = mvals[i]
	
	toplot = melt(round(sqrt(var_mat_m[,1:Nfs,m,p]),3),varnames = c('sigma','f_sigma'))
	toplot$sigma = (rep(sigma_vals,times=(Nfs)))
	toplot$f_sigma =(rep(f_sigma_vals[1:Nfs],each=Ns))
	colnames(toplot)[3] = 'sd'
	toplot$num = as.vector(num_peaks[,1:Nfs,m,p])
	toplot$sd_toshow = rep(sdbreaks_toshow,length.out=Nfs*Ns)
	toplot$num_toshow = rep()
	
contour_sd <- ggplot(toplot, aes(x=f_sigma,y= sigma, z=sd) )+ 
	stat_contour(aes(x=f_sigma,y=sigma,z=sd, colour=as.factor(..level..)),breaks=sdbreaks) + 
	scale_colour_manual(breaks=sdbreaks,values=sdpal) +
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position='none') + 
		 scale_x_continuous(expand=c(0.01,0),limits=range(toplot$f_sigma),breaks=xbreaks,labels=xbreaks)+
		scale_y_continuous(expand=c(0.01,0),limits=ylim,breaks=ybreaks,labels=ybreaks) +xlab('Female stand. dev.') + ylab('Promiscuity')+
		labs(color='Stand. dev.')
		
contour_sd_legend <- ggplot(toplot, aes(x=f_sigma,y= sigma, z=sd_toshow) )+ 
	stat_contour(aes(x=f_sigma,y=sigma,z=sd, colour=as.factor(..level..)),breaks=sdbreaks_toshow) + 
	scale_colour_manual(breaks=sdbreaks_toshow,values=sdpal_toshow) +
	theme_bw() +
	theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position=c(0.6,.53)) + 
	labs(color='Stand. dev.')
		
contour_num <- ggplot(toplot, aes(x=f_sigma,y= sigma, z=num) )+ 
	stat_contour(aes(x=f_sigma,y=sigma,z=num, colour=as.factor(..level..)),breaks=numbreaks) + 
	scale_colour_manual(breaks=numbreaks,values=numpal) +
		theme_bw() +
		theme(text=element_text(family="Helvetica", size=10),plot.title=element_text(size=10) , plot.margin=unit(marg,"cm"),legend.position=c(0.6,.5)) + 
		 scale_x_continuous(expand=c(0.01,0),limits=range(toplot$f_sigma),breaks=xbreaks,labels=xbreaks)+
		scale_y_continuous(expand=c(0.01,0),limits=ylim,breaks=ybreaks,labels=ybreaks) + xlab('Female stand. dev.') + ylab('Promiscuity') + 
		labs(color='# of peaks')

# }

legend_sd <- get_legend(contour_sd_legend)
# contour_sd = contour_sd + theme(legend.position='none') 
legend_num <- get_legend(contour_num)
contour_num = contour_num + theme(legend.position='none') 

contour_tot = arrangeGrob(contour_sd,legend_sd,contour_num,legend_num,nrow=1,widths=c(1,0.3,1,0.3))

	
pdf(file=paste('/Users/eleanorbrush/Documents/research/song_learning_evolution/examples_and_summary_malefath_femalefath_cont_p=',p,'.pdf',sep=''),width=6.83,height=5)
grid.arrange(examples_tot,contour_tot,ncol=1)
dev.off()
