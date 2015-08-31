library(nleqslv)
solve.growth <- function(theta_mat,Dr_mat,R){
	one <- growth_function(theta_mat[,1],Dr_mat[,1],R)
	N_mat <- g_mat <- matrix(ncol=dim(theta_mat)[2],nrow=dim(one)[1])
	N_mat[,1] <- one[,1]
	g_mat[,1] <- one[,2]
	for(i in 2:dim(theta_mat)[2]){
		one <- growth_function(theta_mat[,i],Dr_mat[,i],R)
		N_mat[,i] <- one[,1]
		g_mat[,i] <- one[,2]
	}
	colnames(N_mat) <- colnames(g_mat) <- colnames(theta_mat)
	list(N_mat=N_mat,g_mat=g_mat)
}

plot.growth <- function(growth_mats,filename=NA){
	if(!is.na(filename)){
		png(filename,width=800,height=400)
		par(mfrow=c(1,2))
	}
	xlimit <- sum(growth_mats[[2]][,1]>0)
	for(i in 2:dim(growth_mats[[2]])[2]){
		xlimit <- max(xlimit,max(growth_mats[[1]][which(growth_mats[[2]][,i]>0),i]))
	}
	plot(growth_mats[[2]][,1]~growth_mats[[1]][,1],xlab=paste("N_i(R=",R,")",sep=""),ylab="g_i",ylim=c(0,max(growth_mats[[2]])),xlim=c(min(growth_mats[[1]]),xlimit),col=0)
	for(i in 1:dim(growth_mats[[1]])[2]){
		lines(growth_mats[[1]][,i],growth_mats[[2]][,i],col=1+i)
	}
	plot.new()
	legend("bottomleft",colnames(growth_mats[[1]]),lty=1,col=1:dim(growth_mats[[1]])[2]+1,cex=1)
	if(!is.na(filename)){
		dev.off()
	}
}

optim_exam <- function(R,sp,var="theta",pic=F){##for two species only
	ntest=100
	theta_mat <- result <- Dr_mat <- matrix(rep(0.5,2*ntest),ncol=2)
	theta_mat <- theta_mat*2
	if(var=="theta"){
		var_seq <- theta_mat[,sp] <- seq(0.1,2,length=ntest)
	}else if(var=="Dr"){
		var_seq <- Dr_mat[,sp] <- seq(0.1,0.9,length=ntest)
	}else{
		print("pls provide the correct variable")
	}
	for(i in 1:ntest){
		result[i,] <- growth_function(theta_mat[i,],Dr_mat[i,],R)
	}
	result <- data.frame(result,var_seq)
	if(length(sp)>1){
		sp = "all"
	}
	names(result) <- c("N1_optim","g1_optim",paste(var,sp,sep="_"))
	if(pic){
		plot_optim(result)
	}
	result
}

plot_optim <- function(result,graph=F){
	if(graph){
		filename=paste(names(result)[3],"_optim.png",sep="")
		png(filename,width=1200,height=400)
	}
	scale = 1.8 ## graphing scale
	par(mfrow=c(1,3))
	plot(result[,1]~result[,3],ylab="optimal N1",xlab="",cex.lab=scale,cex.axis=scale,cex=scale)
	plot(result[,2]~result[,3],ylab="optimal g1",xlab=names(result)[3],cex.lab=scale,cex.axis=scale,cex=scale)
	plot(result[,2]/result[,1]~result[,3],ylab="optimal g1/N1",xlab="",cex.lab=scale,cex.axis=scale,cex=scale)
	if(graph){
		dev.off()
	}
}

growth_function <- function(theta_seq,Dr_seq,R){##taking all the other species at equilibrium!
	r_range <- seq(0.1,20,by=0.1)
	i_max1 <- get_ri(theta_seq,Dr_seq,R,r_range)$i_max##narrow down solution span
	r_range <- seq(r_range[min(i_max1)-1],r_range[min(i_max1)+1],by=0.0005)
	solutions <- get_ri(theta_seq,Dr_seq,R,r_range)
	#print(r_range[solutions$i_max])
	if(min(solutions$i_max)==1|max(solutions$i_max)==length(solutions$g_seq)){
		print("touch edge")
	}
	c(solutions$N_seq[solutions$i_max][1],solutions$g_seq[solutions$i_max][1])
}

get_ri <- function(theta_seq,Dr_seq,R,r_range){
	Ng <- c()
	for(i in 1:length(r_range)){
		Ng <- c(Ng,solve.analytical(theta_seq,Dr_seq,R,r_range[i])[1])##only calculate for species 1

	}
	N_seq <- Ng*(2^(Dr_seq[1]*theta_seq[1])*r_range^0.5+r_range+1)/(2^(Dr_seq[1]*theta_seq[1])*r_range^0.5+2*r_range)##only calculate for species 1
	g_seq <- Ng-N_seq
	list(i_max=which(g_seq==max(g_seq)),N_seq=N_seq,g_seq=g_seq)
}

solve.equilibrium <- function(theta_range,Dr_range,theta_controls,Dr_controls,R){
	Dr_mat <- matrix(ncol=length(theta_controls),nrow=length(Dr_range))
	theta_mat <- matrix(ncol=length(Dr_controls),nrow=length(theta_range))
	for(i in 1:length(theta_controls)){
		theta_seq <- rep(theta_controls[i],length(Dr_range))
		Dr_mat[,i] <- solve.analytical(theta_seq,Dr_range,R,1)	
	}
	for(i in 1:length(Dr_controls)){
		Dr_seq <- rep(Dr_controls[i],length(theta_range))
		theta_mat[,i] <- solve.analytical(theta_range,Dr_seq,R,1)	
	}
	rownames(theta_mat) <- theta_range
	rownames(Dr_mat) <- Dr_range
	colnames(Dr_mat) <- paste("theta=",theta_controls,sep="")
	colnames(theta_mat) <- paste("Dr=",Dr_controls,sep="")
	list(theta_mat=theta_mat,Dr_mat=Dr_mat)
}

plot.equilibrium <- function(equi_mats,filename=NA){
	if(!is.na(filename)){
		png(filename,width=800,height=400)
		par(mfrow=c(1,2))
	}
	names <- c("theta","Dr")
	position <- c("topright","topleft")
	for(i in 1:2){
		plot(equi_mats[[i]][,1]~as.numeric(rownames(equi_mats[[i]])),ylim=range(equi_mats[[i]]),xlab=paste(names[i],"(R=",R,")",sep=""),ylab="Ne",col=0)
		for(j in 1:dim(equi_mats[[i]])[2]){
			lines(as.numeric(rownames(equi_mats[[i]])),equi_mats[[i]][,j],col=j+1)
		}
		legend(position[i],colnames(equi_mats[[i]]),lty=1,col=1:dim(equi_mats[[i]])[2]+1)
	}
	if(!is.na(filename)){
		dev.off()
	}
}

solve.analytical <- function(theta_seq,Dr_seq,R,r){##given r solve for N+g; r=1 is the special case for equilibrium where N+g=Ne
	r_seq <- c(r,rep(1,length(theta_seq)-length(r)))
	C <- get.C(theta_seq,Dr_seq,R,r_seq)
	analytical.sol(C,theta_seq,Dr_seq,r_seq)
}

get.C <- function(theta_seq,Dr_seq,R,r_seq){
	init.C <- 1/mean(theta_seq)
	C <- nleqslv(init.C,R_constraint,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R,r_seq=r_seq)
	if(C[[3]]!=1&&C[[3]]!=2){
		print(C[[4]])
	}
	C[[1]]
}

R_constraint <- function(C,theta_seq,Dr_seq,R,r_seq){
	R-sum(theta_seq*analytical.sol(C,theta_seq,Dr_seq,r_seq))
}

analytical.sol <- function(C,theta_seq,Dr_seq,r_seq){
	2/exp(1)*(C*theta_seq*r_seq^(0.5/theta_seq))^(1/(Dr_seq-1))
}

SAD_shape <- function(Sp,R,controls,vary="theta",distr="normal",scale=1){
	##first create the series of variable values ranging from 0 to 2 for different distribution types
	sample <- seq(0,2,length=Sp+2)[-c(1,Sp+2)]
	if(distr=="normal"){
		density <- dnorm(sample,1,0.3)##because pnorm(0,1,0.3)<0.001 & pnorm(2,1,0.3)>0.999
	}else if(distr=="exponential"){
		density <- dexp(sample,4)##because pexp(0,4)<0.001 & pexp(2,4)>0.999
	}else if(distr=="uniform"){
		density <- 0.5
	}
	density <- density/sum(density) ##rescale to sum to 1
	Ne_mat <- pdf_mat <- matrix(nrow=Sp,ncol=length(controls))
	for(i in 1:length(controls)){
		if(vary=="theta"){
			control="Dr"
			Ne_mat[,i] <- sort(solve.analytical(sample,rep(controls[i],Sp),R,1))
			pdf_mat[,i] <- abs(density*var_N_de(Ne_mat[,i],sample,controls[i],vary))
			pdf_mat[,i] <- pdf_mat[,i]/sum(pdf_mat[,i]) ##rescale to sum to 1
		}else if(vary=="Dr"){
			control="theta"
			Ne_mat[,i] <- sort(solve.analytical(rep(controls[i],Sp),sample/2.5,R,1))## assumed to be half the normal distribution of theta
			pdf_mat[,i] <- abs(density*var_N_de(Ne_mat[,i],sample/2.5,controls[i],vary))
			pdf_mat[,i] <- pdf_mat[,i]/sum(pdf_mat[,i]) 
		}
	}
		plot(pdf_mat[,1]~Ne_mat[,1],xlim=c(min(Ne_mat),max(Ne_mat[,1])/5),ylim=range(pdf_mat),xlab=paste("abundance(",vary,distr,")"),ylab="pdf",col=0,cex=scale,cex.lab=scale,cex.axis=scale)
	for(i in 1:length(controls)){
		lines(Ne_mat[,i],pdf_mat[,i],col=i+1,lwd=scale)
	}	
}

var_N_de <- function(Ne_seq,var_seq,control,vary="theta"){
	if(vary=="theta"){
		var_seq*(control-1)/Ne_seq
	}else if(vary=="Dr"){
		(var_seq-1)/log(exp(2)*Ne_seq/2)/Ne_seq
	}else{
		print("Please provide the correct variable")
	}
}