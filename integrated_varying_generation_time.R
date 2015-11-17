library(nleqslv)

multi_predator_1prey <- function(inits,thetas,Drs,rs,Gs,R0,nstep){
	Nseq <- matrix(nrow=nstep,ncol=length(thetas))
	Nseq[1,] <- inits
	for(t in 2:nstep){
		clocks <- cumu_growths <- rep(0,length=length(thetas))
		cumu_growths[-1] <- cumu_growths[-1] + multi_growth(Nseq[t-1,-1],thetas[-1],Drs[-1],rs[-1],Nseq[t-1,1]*thetas[1])
		cumu_growths[1] <- multi_growth(Nseq[t-1,1],thetas[1],Drs[1],rs[1],R0)- sum(Nseq[t-1,-1]*thetas[-1])/thetas[1]#lag in predation
		clocks <- clocks+1
		for(sp in 1:length(thetas)){
			if(clocks[sp]<Gs[sp]){#during the current generation, abundance doesn't change
				Nseq[t,sp] <- Nseq[t-1,sp]
			}else{#when the current generation ends, update abundance and reset clock and cumulated growth
				Nseq[t,sp] <- Nseq[t-1,sp]+cumu_growths[sp]
				clocks[sp] <- 0
				cumu_growths[sp] <- 0
			}
		}		
	}

	time <- 1: nstep
	plot(Nseq[,1]~time,col=1,ylim=range(Nseq),ylab="abundance")
	lines(time,Nseq[,1],col=1)
	for(i in 2:dim(Nseq)[2]){
		points(time,Nseq[,i],col=i)
		lines(time,Nseq[,i],col=i)
	}
	#legend("bottomright",c("prey",paste("predator",1:(dim(Nseq)[2]-1))),col=1:length(thetas),pch=1,lty=1)
	Nseq
}

multi_integrated <- function(init.Ns0,thetas,Drs,rs,Gs,R0,nstep){
	Nseq <- matrix(nrow=nstep,ncol=length(thetas))
	Nseq[1,] <- init.Ns0
	clocks <- cumu_growths <- rep(0,length=length(thetas))
	for(i in 2:nstep){
		cumu_growths <- cumu_growths + multi_growth(Nseq[i-1,],thetas,Drs,rs,R0)
		clocks <- clocks+1
		for(sp in 1:length(thetas)){
			if(clocks[sp]<Gs[sp]){#during the current generation, abundance doesn't change
				Nseq[i,sp] <- Nseq[i-1,sp]
			}else{#when the current generation ends, update abundance and reset clock and cumulated growth
				Nseq[i,sp] <- Nseq[i-1,sp]+cumu_growths[sp]
				clocks[sp] <- 0
				cumu_growths[sp] <- 0
			}
		}
	}
	time <- 1:nstep
	for(sp in 1:length(Drs)){
		if(sp==1){
			plot(Nseq[,sp]~time,col=sp+1,ylim=range(Nseq),ylab="abundance",xlab=paste("time (Dr=",mean(Drs),")"))
		}else{
			points(time,Nseq[,sp],col=sp+1)
		}
		lines(time,Nseq[,sp],col=sp+1)
	}
	Nseq
}

multi_growth <- function(Ns,thetas,Drs,rs,R0){
	Cs <- (rs+1)*Ns*thetas
	NDs <- Ns^Drs
	Rs.init <- NDs/sum(NDs)*R0*sum(Cs)/(R0+sum(Cs))
	Rs <- nleqslv(Rs.init,dynamic_constraints,NDs=NDs,Cs=Cs,R0=R0)
	if(Rs[[3]]!=1){
		print(Rs[[4]])
	}
	Rs <- Rs[[1]]
	Rs/thetas-Ns
}

dynamic_constraints <- function(Rs,NDs,Cs,R0){
	Rs^2 - (Cs-Rs)*(R0-sum(Rs))*NDs
}