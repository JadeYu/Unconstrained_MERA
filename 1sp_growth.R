library(nleqslv)
single_growth <- function(r,R0,N,Dr,theta){
	if(N<R0/r){
		r*N^Dr/theta/(N^Dr+1)-1
	}else{
		R0*N^(Dr-1)/theta/(N^Dr+1)-1
	}
}

single_growth <- function(r,R0,N,Dr,theta){##the stock version
	r*(R0-theta*N)*(r*N)^(Dr-1)/(theta*((r*N)^Dr+1))
}

SSN_eq <- function(init,R0,Dr,theta){
	R0*init^(Dr-1)-theta*(init^Dr+1)
}

solve_SSN <- function(R0,Dr,theta){
	init <- mean(c(R0/theta/2,R0/theta-1))
	SSN <- nleqslv(init,SSN_eq,R0=R0,Dr=Dr,theta=theta)
	if(SSN[[3]]!=1){
		print(SSN[[4]])
	}
	SSN[[1]]
}

plot_growth <- function(SSN,r,R0,Dr,theta,ng=1,comp=T){
	N_seq <- seq(10,R0,length=50)
	g_seq <- unlist(lapply(N_seq,single_growth,r=r,R0=R0,Dr=Dr,theta=theta))
	if(ng==1){
		plot(g_seq~N_seq,xlab="abundance N",ylab="net growth rate g",col=2,pch=16)
	}else{
		points(N_seq,g_seq,pch=16,col=ng+1)
	}
	if(comp){
		logi_g <- (r/2/theta-1)*(1-N_seq/SSN)
		points(N_seq,logi_g,pch=16)
	}
}

plot_growth <- function(r,R0,Dr,theta,ng=1,comp=T){
	SSN <- R0/theta
	N_seq <- seq(1,R0,length=50)
	g_seq <- unlist(lapply(N_seq,single_growth,r=r,R0=R0,Dr=Dr,theta=theta))
	if(ng==1){
		plot(g_seq~N_seq,xlab="abundance N",ylab="net growth rate g",col=2,pch=16)
	}else{
		points(N_seq,g_seq,pch=16,col=ng+1)
	}
	if(comp){
		logi_g <- r*(1-N_seq/SSN)
		points(N_seq,logi_g,pch=16)
	}
}

plot_dynamic <- function(SSN,r,R0,Dr,theta,nstep,ng=1,comp=T){
	N <- 1
	N_seq <- N
	t <- 1
	while (t <nstep){
		N <- N*(1+single_growth(r,R0,N,Dr,theta))
		N_seq <- c(N_seq,N)
		t <- t+1
	}
	time <- (1:nstep)
	if(ng==1){
		plot(N_seq~time,xlab="time",ylab="abundance N",col=2,pch=16)
	}else{
		points(time,N_seq,pch=16,col=ng+1)
	}
	if(comp){
		t <- 1
		N <- 1
		N_comp <- N
		while (t < nstep){
			N <- N*(1+(r/2/theta-1)*(1-N/SSN))
			N_comp <- c(N_comp,N)
			t <- t+1
		}
		points(time,N_comp,pch=16)
	}
}

growth_function <- function(r,R0,Dr_seq,theta,nstep){
	par(mfrow=c(1,2))
	for(i in 1:length(Dr_seq)){
		SSN <- solve_SSN(R0,Dr_seq[i],theta)
		if(i ==1){
			plot_growth(SSN,r,R0,Dr_seq[i],theta)
		}else{
			plot_growth(SSN,r,R0,Dr_seq[i],theta,ng=i,comp=F)
		}	
	}
	for(i in 1:length(Dr_seq)){
		SSN <- solve_SSN(R0,Dr_seq[i],theta)
		if(i ==1){
			plot_dynamic(SSN,r,R0,Dr_seq[i],theta,nstep,comp=F)
		}else{
			plot_dynamic(SSN,r,R0,Dr_seq[i],theta,nstep,ng=i,comp=F)
		}	
	}
	
	legend("bottomright",c(paste("Dr=",Dr_seq),"logistic growth"),pch=16,col=c(1:length(Dr_seq)+1,1))
}
