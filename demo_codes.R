MERA_growth <- function(N_seq,theta_seq,Dr_seq,R){
	C <- sum(N_seq^Dr_seq)
	N_seq^Dr_seq*R/theta_seq/(1+C)-N_seq
}

growth_1sp <- function(init.Ns,r,theta_seq,Dr_seq,R,nstep){
	SSN <- get_SSN(theta_seq,Dr_seq,R)
	N_seq <- init.Ns
	N_mat <- matrix(nrow=nstep,ncol=length(init.Ns))
	for(t in 1:nstep){
		N_mat[t,] <- N_seq
		potentials <- r*N_seq
		G <- MERA_growth(N_seq,theta_seq,Dr_seq,R)
		for(i in 1:length(init.Ns)){
			N_seq[i] <- N_seq[i]+ min(potentials[i],G[i])
		}
	}
	plot(N_mat[,1],xlab="time",ylab="abundance",ylim=range(N_mat))
	if(length(init.Ns)>1){
		for(i in 2:length(init.Ns)){
		points(N_mat[,i],col=i)
		}	
	}
	N_mat
}

growth_1sp <- function(init.Ns,r=2,theta_seq,Dr_seq,R,nstep){
	SSN <- get_SSN(theta_seq,Dr_seq,R)
	N_seq <- init.Ns
	N_mat <- matrix(nrow=nstep,ncol=length(init.Ns))
	for(t in 1:nstep){
		N_mat[t,] <- N_seq
		R_avail <- min(r*sum(N_seq)^2,R)## total resource available grows linearly to the total resource already obtained
		G <- MERA_growth(N_seq,theta_seq,Dr_seq,R_avail)
		N_seq <- N_seq+ G
	}
	plot(N_mat[,1],xlab="time",ylab="abundance",ylim=range(N_mat))
	if(length(init.Ns)>1){
		for(i in 2:length(init.Ns)){
		points(N_mat[,i],col=i)
		}	
	}
	N_mat
}

library(nleqslv)
get_SSN <-function(theta_seq,Dr_seq,R){
	inits <- R/theta_seq/2
	SSN <- nleqslv(inits,MERA_growth,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R)
	if(SSN[[3]]!=1){
		print(SSN[[4]])
	}
	SSN[[1]]
}

G <- MERA_growth(N_seq,theta_seq,Dr_seq,R)



growth_1sp(12,r=0.15,theta_seq=1,Dr_seq=0.3,R=100,nstep=30)

theta_seq <- c(1,1,1)
Dr_seq <- c(0.1,0.1,0.5)
R <- 100

SSN <- get_SSN(theta_seq,Dr_seq,R)

growth_1sp(init.Ns=c(5,1,10),r=1.5,theta_seq,Dr_seq,R=100,nstep=50)