library(nleqslv)
init.Ns <- c(5,10,20)
rs <- c(100,100,100)
thetas <- c(3,2,1)
Drs <- c(0.9,0.5,0.1)
R0 <- 100
nstep <- 10

init.Ns <- round(rnorm(10,100,10))
rs <- rep(100,10)
thetas <- rnorm(10,5,2)
Drs <- runif(10,0.7,1)
R0 <- 10000
nstep <- 100

multi_competitor(init.Ns,rs,thetas,Drs,R0,nstep)

multi_competitor <- function(init.Ns,rs,thetas,Drs,R0,nstep){
	Nseq <- matrix(ncol=length(Drs),nrow=nstep)
	Nseq[1,] <- init.Ns
	for(t in 2:nstep){
		Nseq[t,] <- update(Nseq[t-1,],rs,thetas,Drs,R0)
	}
	time <- 1:nstep
	for(sp in 1:length(Drs)){
		if(sp==1){
			plot(Nseq[,sp]~time,col=sp+1,ylim=range(Nseq))
		}else{
			points(time,Nseq[,sp],col=sp+1)
		}
	}
	Nseq
}

update <- function(init.Ns,rs,thetas,Drs,R0){
	Rm <- sum((rs+1)*init.Ns*thetas)
	Rr <- Rm*R0/(Rm+R0)
	gs <- MERA.growth(init.Ns,thetas,Drs,Rr)
	init.Ns*(1+gs)
}

MERA.growth <- function(Ns,thetas,Drs,R0){## growth function 
	PR <- Ns^(Drs)
	(R0 * PR/sum(PR)/thetas - Ns)/Ns
}

MERA.SSN <- function(thetas,Drs,R0){
	C <- get.C(thetas,Drs,R0) 
	(C*thetas)^(1/(Drs-1))
}

get.C <- function(thetas,Drs,R0){
	init.C <- (R0/sum(thetas^(Drs/(Drs-1))))^(mean(Drs)-1)
	C <- nleqslv(init.C,constraint,thetas=thetas,Drs=Drs,R0=R0)
	if(C[[3]]!=1){
		print(C[[4]])
	}
	C[[1]]
} 

constraint <- function(C,thetas,Drs,R0){
	R0 - sum(C^(1/(Drs-1))*thetas^(Drs/(Drs-1)))
}