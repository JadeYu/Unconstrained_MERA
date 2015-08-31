analytical.sol <- function(C,theta_seq,Dr_seq,r_seq){
	2/exp(1)*(C*theta_seq*r_seq^(0.5/theta_seq))^(1/(Dr_seq-1))
}

R_constraint <- function(C,theta_seq,Dr_seq,R,r_seq){
	R-sum(theta_seq*analytical.sol(C,theta_seq,Dr_seq,r_seq))
}

solving_growth <- function(x,N_seq,theta_seq,Dr_seq,R){
	r_seq <- x[1:(length(x)-1)]
	C <- x[length(x)]
	y <- numeric(length(x))
	Nseq <- analytical.sol(C,theta_seq,Dr_seq,r_seq)*(2^(Dr_seq*theta_seq)*r_seq^0.5+r_seq+1)/(2^(Dr_seq*theta_seq)*r_seq^0.5+2*r_seq)
	y[1:(length(x)-1)] <- N_seq- Nseq
	y[length(x)] <- R_constraint(C,theta_seq,Dr_seq,R,r_seq)
	y
}

solving_growth_Dr1 <- function(C,N_seq,theta_seq,R){
	Ng <- growth_Dr1(C,N_seq,theta_seq)
	R-sum(theta_seq*Ng)
}##not working yet

growth_Dr1 <- function(C,N_seq,theta_seq){
	N_seq*(2*C^(2*theta_seq)+factorial(2*theta_seq)/factorial(theta_seq)*C^theta_seq)/(C^(2*theta_seq)+factorial(2*theta_seq)/factorial(theta_seq)*C^theta_seq+factorial(2*theta_seq))
}##not working yet


numeric.growth <- function(N_seq,theta_seq,Dr_seq,R,var="N"){
	if(sum(Dr_seq==1)){
		init.C <- get.C(theta_seq,rep(0.95,length(N_seq)),R,rep(1,length(N_seq)))
		sol <- nleqslv(init.C,solving_growth_Dr1,N_seq=N_seq,theta_seq=theta_seq,R=R,control=list(btol=10^-6,maxit=1000))
		if(sol[[3]]!=1){
			print(sol[[4]])
			return()
		}
		C <- sol[[1]]
		N_next <- growth_Dr1(C,N_seq,theta_seq)
	}else{
		init.r <- rep(1,length(N_seq))
		init.C <- get.C(theta_seq,Dr_seq,R,init.r)
		init.x <- c(init.r,init.C)
		solution <- nleqslv(init.x,solving_growth,N_seq=N_seq,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R,control=list(xtol=10^-15))
		if(solution[[3]]!=1){
			print(solution[[4]])
			return()
		}
		C <- solution[[1]][(length(N_seq)+1)]
		r_seq <- solution[[1]][1:length(N_seq)]
		N_next <- analytical.sol(C,theta_seq,Dr_seq,r_seq)
	}
	if(var=="N"){
		result <- N_next
	}else{
		result <- (N_next-N_seq)/N_seq
	}
	result
}

dynamics <- function(N_seq,R,Dr_seq,theta_seq,nstep){
	time_track <- matrix(nrow=nstep,ncol=length(N_seq))
	time_track[1,] <- N_seq
	for(i in 2:nstep){
		time_track[i,] <- numeric.growth(time_track[i-1,],theta_seq/10,Dr_seq,R/10)
	}
	colnames(time_track) <- paste("sp",1:length(N_seq),sep="")
	time_track
}

get.C <- function(theta_seq,Dr_seq,R,r_seq){
	init.C <- 1/mean(theta_seq)
	C <- nleqslv(init.C,R_constraint,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R,r_seq=r_seq)
	if(C[[3]]!=1&&C[[3]]!=2){
		print(C[[4]])
	}
	C[[1]]
}