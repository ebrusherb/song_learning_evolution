
sigma2 = 0
sigmay2_init = 1.8
sigmax2_init = 1
rho_init = 0

analytical = array(0,dim=c(9,4))
analytical[3,1]=sigmay2_init-sigma2
analytical[3,2]=sigmay2_init
analytical[4,1]=sigmax2_init
analytical[5,1]=sigmax2_init
analytical[5,2]=sigmax2_init*(sigma2+sigmax2_init)/(2*sigmax2_init+sigma2)
analytical[6,1]=sigmax2_init
analytical[6,2]=sigmay2_init
sigmax2_star = (3*sigmay2_init-5*sigma2+sqrt(sigma2^2+sigmay2_init*(9*sigmay2_init-30*sigma2)))/6
analytical[9,1]=sigmax2_star
analytical[9,2]=sigmay2_init
analytical[9,3]=sigmax2_star*sigmay2_init/(sigma2+sigmax2_star)

# #########

# # sigma2 = 0
# sigmay2_init = 1.8
# sigmax2_init = 1
# rho_init = 0

# steps = 30

# variables_mat <- recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)

# my_cols = rainbow(9)
# layout(matrix(1:3,ncol=3))
# plot(1:(steps+1),1:(steps+1),type='n',ylim=c(0,sigmay2_init))
# for(i in 1:9){
	# points(1:(steps+1),variables_mat[i,1,],col=my_cols[i],t='o')
# }
# legend(x=(0.5)*steps,y=1.5,legend=1:9,pch=1,col=my_cols)
# plot(1:(steps+1),1:(steps+1),type='n',ylim=c(0,sigmay2_init))
# for(i in 1:9){
	# points(1:(steps+1),variables_mat[i,2,],col=my_cols[i],t='o')
# }
# plot(1:(steps+1),1:(steps+1),type='n',ylim=c(0,1))
# for(i in 1:9){
	# points(1:(steps+1),variables_mat[i,4,],col=my_cols[i],t='o')
# }


#######
sigma2_vals = seq(0,10,length.out=50)
sigmay2_init_vals = seq(0,10,length.out=25)
sigmax2_init = 2
rho_init = 0.7
steps = 1000
transient = 20

x_sigma2 = length(sigma2_vals)
x_sigmay2 = length(sigmay2_init_vals)

variables_transient = array(NA,dim=c(9,5,x_sigma2,x_sigmay2))
variables_final = array(NA,dim=c(9,5,x_sigma2,x_sigmay2))

for(i in 1:x_sigma2){
	for(j in 1:x_sigmay2){
		sigma2 = sigma2_vals[i]
		sigmay2_init = sigmay2_init_vals[j]
		variables_mat <- recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
		variables_transient[,,i,j]=variables_mat[,,transient]
		variables_final[,,i,j]= variables_mat[,,steps+1]
	}
}

inf_cutoff = 1000
variables_toplot = variables_transient
variables_toplot[which(variables_final>=inf_cutoff)]=inf_cutoff

# variables_toplot = variables_final

not_boring = c(1:3,4:6,7:9)
k = 1

breaks=c(0,10^seq(-5,2,length.out=20),max(variables_toplot[,k,,][which(!is.nan(variables_toplot[,k,,]))]))
layout(matrix(1:length(not_boring),nrow=3))
for(i in not_boring){
	image(sigma2_vals,sigmay2_init_vals,variables_toplot[i,k,,],breaks=breaks,col=c(heat.colors(length(breaks)-2),'black'))
	if(i==3){
		abline(0,1)
	}
	if(i==9){
		abline(0,(30+sqrt(864))/18)
	}
}

# mycols = rainbow(9)
# layout(matrix(1:(x_sigmay2),ncol=3,byrow=TRUE))
# for(j in 1:x_sigmay2){
	# ymax = min(max((variables_transient[not_boring,k,,j])[which(!is.nan(variables_transient[not_boring,k,,j]))]),2)
	# plot(sigma2_vals,1:x_sigma2,type='n',ylim=c(0,ymax))
	# for(i in not_boring){
		# points(sigma2_vals,variables_transient[i,k,,j],col=my_cols[i],t='o')
	# }
	# abline(h=sigmay2_init_vals[j])
	# abline(v=sigmay2_init_vals[j])
	# if(j ==1){
		# legend(x=mean(sigma2_vals),y=ymax,legend=(1:9)[not_boring],pch=1,col=my_cols[not_boring])}
# }

# i=3
# j=9
# sigma2 = sigma2_vals[i]
# sigmay2_init = sigmay2_init_vals[j]
# sigmax2_star = (3*sigmay2_init-5*sigma2_vals+sqrt(sigma2_vals^2+sigmay2_init*(9*sigmay2_init-30*sigma2_vals)))/6