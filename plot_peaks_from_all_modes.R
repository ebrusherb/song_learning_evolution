par(mfrow=c(2,4))

t0 = 1
tseq = seq(100,1000,by=200)

ylim = c(0,0.15)
xlim = c(73,133)

plot(d3$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d3$Pm[,t])
}

lines(d3$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d3$Pf[,t],col='red')
}

plot(apply(d4$Pm[,,t0],1,sum),t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(apply(d4$Pm[,,t],1,sum))
}

lines(d4$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d4$Pf[,t],col='red')
}

plot(apply(d6$Pm[,,t0],1,sum),t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(apply(d6$Pm[,,t],1,sum))
}

lines(d6$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d6$Pf[,t],col='red')
}

plot(d7$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d7$Pm[,t])
}

lines(d7$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d7$Pf[,t],col='red')
}

plot(d8$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d8$Pm[,t])
}

lines(apply(d8$Pf[,,t0],2,sum),t='l',col='red')
for(t in tseq){
	lines(apply(d8$Pf[,,t],2,sum),col='red')
}

plot(d9$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d9$Pm[,t])
}

lines(d9$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d9$Pf[,t],col='red')
}

plot(d11$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d11$Pm[,t])
}

lines(apply(d11$Pf[,,t0],2,sum),t='l',col='red')
for(t in tseq){
	lines(apply(d11$Pf[,,t],2,sum),col='red')
}

plot(d12$Pm[,t0],t='l',ylim=ylim,xlim=xlim)
for(t in tseq){
	lines(d12$Pm[,t])
}

lines(d12$Pf[,t0],t='l',col='red')
for(t in tseq){
	lines(d12$Pf[,t],col='red')
}



