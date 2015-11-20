library(nleqslv)

test_continuity <- function(init.Ns0,I_thetas,Drs,I_rs,Gs_abs,I_R0,dt_seq){
	for(i in 1:length(dt_seq)){
		Nseq <- cross_time(init.Ns0,I_thetas,Drs,I_rs,Gs=Gs_abs/dt_seq[i],I_R0,round(50/dt_seq[i]),dt_seq[i],F)
		plot_dynamics(Nseq,dt,"dt",dt_seq[i],pch=i,new=(i==1))
	}
}

cross_time <- function(init.Ns0,I_thetas,Drs,I_rs,Gs,I_R0,nstep,dt,graph){
	R0 <- I_R0*dt
	thetas <- I_thetas*dt
	rs <- exp(I_rs*dt)-1 ##transfer into discrete growth rates
	Nmat <- multi_integrated(init.Ns0,thetas,Drs,rs,Gs,R0,nstep,dt,graph)
	Nmat 
}

multi_integrated <- function(init.Ns0,thetas,Drs,rs,Gs,R0,nstep,dt,graph=T){
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
	if(graph){
		plot_dynamics(Nseq,dt,"r",mean(rs),1)
	}
	Nseq
}

plot_dynamics <- function(Nseq,dt,xlab,values,pch,new=T){
	time <- 1:dim(Nseq)[1]*dt
	for(sp in 1:dim(Nseq)[2]){
		if(sp==1&new){
			plot(Nseq[,sp]~time,col=sp+1,ylim=range(Nseq),ylab="abundance",xlab=paste("time (",xlab,"=",values,")",sep=""),pch=pch)
		}else{
			points(time,Nseq[,sp],col=sp+1,pch=pch)
		}
		lines(time,Nseq[,sp],col=sp+1)
	}

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