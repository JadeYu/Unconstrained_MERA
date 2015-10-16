growth <- function(r,N,SSN,partial=T){
	if(partial){
		r*(1-N/SSN)/(1+(r-1)*N/SSN)
	}else{
		r*(1-N/SSN)/(1+r*N/SSN)
	}
}

logi_growth <- function(r,N,SSN){
	r*(1-N/SSN)
}

plot_dynamic <- function(r,SSN,nstep,partial=T){
	N1_seq <- N2_seq <- 1
	t <- 1
	while (t < nstep){
		N1_seq <- c(N1_seq,N1_seq[length(N1_seq)]*(1+growth(r,N1_seq[length(N1_seq)],SSN,partial)))
		N2_seq <- c(N2_seq,N2_seq[length(N2_seq)]*(1+logi_growth(r,N2_seq[length(N2_seq)],SSN)))
		t <- t+1
	}
	time <- (1:nstep)
	plot(N2_seq~time,xlab=paste("r=",r,sep=""),ylab="abundance N",col=2,pch=16)
	points(time,N1_seq,pch=16,col=1)
}
par(mfrow=c(2,3))
SSN <- 100
r <- 1.3
nstep=20

r_seq <- c(0.5,1,1.5,2,3)
ns_seq <- c(30,30,20,20,20)

for(i in 1:5){
	plot_dynamic(r_seq[i],SSN,ns_seq[i])
}

for(i in 1:5){
	plot_dynamic(r_seq[i],SSN,ns_seq[i],partial=F)
}

plot.new()
legend("bottomleft",c("MERA","logistic"),pch=16,col=c(1,2))
