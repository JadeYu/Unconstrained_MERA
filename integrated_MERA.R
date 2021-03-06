init.Ns <- c(5,50,100)
rs <- c(0.08,0.05,0.02)
thetas <- c(1,1,1)
Drs <- c(0.1,0.1,0.1)
Gs_abs <- c(40,30,30)
R0 <- 10000
nstep <- 1000

multi_integrated(init.Ns,thetas,Drs,rs,Gs,R0,nstep,T)

#transform into instant values
I_thetas=thetas
I_rs=log(1+rs)
I_R0=R0
dt = 0.5

cross_time(init.Ns0=init.Ns,I_thetas=I_thetas,Drs=Drs,I_rs=I_rs,Gs=Gs_abs/dt,I_R0=I_R0,nstep=nstep,dt=dt,graph=T)


plot.new()
legend("bottomleft",paste("species",1:3),col=2:4,pch=1)

##single species
#test continuity
init.Ns=5
Drs=0

#parameter values when dt=1
thetas=1
<<<<<<< HEAD
rs_logi = 15
=======
rs = 1
R0=270000
>>>>>>> 4998aa0699fd71d3d3c9a004d7cdec3a6d9b084a
#transform into instant values
I_thetas=thetas
Gs_abs=1
I_rs=log(1+rs_logi)/Gs_abs
SSN = 100


#time frame and step length
<<<<<<< HEAD
nstep = 50000
dt = 0.001
rs = exp(I_rs*dt)-1
I_R0 = I_thetas*(rs+1)/rs*SSN

cross_time(init.Ns0=init.Ns,I_thetas=I_thetas,Drs=Drs,I_rs=I_rs,Gs=Gs_abs/dt,I_R0=I_R0,nstep=nstep,dt=dt,graph=T)

logi_plot(SSN,nstep,init.Ns,rs_logi,Gs_abs)
=======
nstep = 100
dt = 0.2

Gs_abs = 1
dt_seq = c(0.1,0.5,1)
test_continuity(init.Ns,I_thetas,Drs,I_rs,Gs_abs,I_R0,dt_seq)

cross_time(init.Ns0=init.Ns,I_thetas=I_thetas,Drs=Drs,I_rs=I_rs,Gs=Gs_abs,I_R0=I_R0,nstep=nstep/dt,dt=dt,graph=T)
>>>>>>> 4998aa0699fd71d3d3c9a004d7cdec3a6d9b084a


multi_integrated(init.Ns=5,thetas=1,Drs=0,rs=0.4,Gs=1,R0=500,nstep=500,T)
logi_growth(5,r=1,SSN=90.48,nstep=20)

plot.new()
legend("bottomright",c("MERA","logistic"),pch=1,col=2:1)

init.Ns <- c(5,50,100)
init.Ns1 <- c(1,45,50,105)
rs <- c(0.5,0.3,0.1)
thetas <- c(1,1,1)
Drs <- c(0.5,0.5,0.5)
R0 <- 200
nstep <- 50
multi_integrated(init.Ns,thetas,Drs,rs,R0,nstep)
legend("bottomleft",paste("species",1:3),col=2:4,pch=1)

get_SSN(thetas,Drs,rs,R0)

##multiple predators
inits <- c(1000,5,50,100)
rs <- c(10,0.5,0.3,0.1)
thetas <- c(1,1,1,1)
Drs <- c(0,0.5,0.5,0.5)
R0 <- 5000
nstep <- 50

multi_predator_1prey(inits,thetas,Drs,rs,R0,nstep)

plot.new()
legend("bottomleft",c("prey",paste("predator",1:(length(thetas)-1))),col=1:length(thetas),pch=1,lty=1)

library(nleqslv)



multi_integrated <- function(init.Ns0,thetas,Drs,rs,R0,nstep){
	Nseq <- matrix(nrow=nstep,ncol=length(thetas))
	Nseq[1,] <- init.Ns0
	for(i in 2:nstep){
		Nseq[i,] <- multi_growth(Nseq[i-1,],thetas,Drs,rs,R0)
	}
	time <- 1:nstep
	for(sp in 1:length(Drs)){
		if(sp==1){
			plot(Nseq[,sp]~time,col=sp+1,ylim=range(Nseq),ylab="abundance",xlab=paste("time (Dr=",Drs[1],")"))
		}else{
			points(time,Nseq[,sp],col=sp+1)
		}
		lines(time,Nseq[,sp],col=sp+1)
	}
	Nseq
}

multi_predator_1prey <- function(inits,thetas,Drs,rs,R0,nstep){
	Nseq <- matrix(nrow=nstep,ncol=length(thetas))
	Nseq[1,] <- inits
	for(t in 2:nstep){
		Nseq[t,-1] <- multi_growth(Nseq[t-1,-1],thetas[-1],Drs[-1],rs[-1],Nseq[t-1,1]*thetas[1])
		Nseq[t,1] <- multi_growth(Nseq[t-1,1],thetas[1],Drs[1],rs[1],R0)- sum(Nseq[t-1,-1]*thetas[-1])/thetas[1]	 #time lag
		#Nseq[t,1] <- multi_growth(Nseq[t-1,1],thetas[1],Drs[1],rs[1],R0)- sum(Nseq[t,-1]*thetas[-1])/thetas[1]	#No time lag
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

multi_growth <- function(Ns,thetas,Drs,rs,R0){
	Cs <- (rs+1)*Ns*thetas
	NDs <- Ns^Drs
	Rs.init <- NDs/sum(NDs)*R0*sum(Cs)/(R0+sum(Cs))
	Rs <- nleqslv(Rs.init,dynamic_constraints,NDs=NDs,Cs=Cs,R0=R0)
	if(Rs[[3]]!=1){
		print(Rs[[4]])
	}
	Rs <- Rs[[1]]
	Rs/thetas
}

dynamic_constraints <- function(Rs,NDs,Cs,R0){
	Rs^2 - (Cs-Rs)*(R0-sum(Rs))*NDs
}

get_SSN <- function(thetas,Drs,rs,R0){
	Ru.init <- (sum((thetas/rs)^(1/(Drs-1))*thetas)/R0)^(mean(Drs)-1)
	Ru <- nleqslv(Ru.init,SS_constraint,thetas=thetas,Drs=Drs,rs=rs,R0=R0)
	if(Ru[[3]]!=1){
		print(Ru[[4]])
	}
	Ru <- Ru[[1]]
	(thetas/rs/Ru)^(1/(Drs-1))
}

SS_constraint <- function(Ru,thetas,Drs,rs,R0){
	R0-Ru-sum((thetas/rs/Ru)^(1/(Drs-1))*thetas)
}

logi_growth <- function(init,r,SSN,nstep){
	Nseq <- numeric(nstep)
	Nseq[1] = init
	for(i in 2:nstep){
		Nseq[i] = Nseq[i-1]*(1+ r* (1-Nseq[i-1]/SSN))
	}
	time <- 1:nstep
	points(time,Nseq)
}