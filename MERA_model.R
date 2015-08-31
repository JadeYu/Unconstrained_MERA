######I. MERA model #########################################
library(nleqslv)
##functions
de_series_ul_cp<- function(x,C_seq,R,dr_seq,re_seq){
	b_seq <- x[1:(length(x)/2)]
	d_seq <- x[((length(x)/2)+1):length(x)]
	Cprime_seq <- C_seq+b_seq-d_seq
	R_resi <- R-crossprod(Cprime_seq,re_seq)
	y <- numeric()
	for(i in 1:length(C_seq)){
		y[(2*i-1):(2*i)] <- DE_ul(C_seq[i],re_seq[i],b_seq[i],d_seq[i],dr_seq[i],R_resi)
	}
	y
}

DE_ul <- function(C,re,b,d,dr,R_resi){
	y <- numeric()
	y[1] <- (4*re/exp(1))^(re*dr)*b*(re*(C+b-d))^(re*(1-dr))-(C-b-d)*R_resi^re
	y[2] <- d*R_resi^re-(re/exp(1))^(re*dr)*(C-b-d)*(re*(C+b-d))^(re*(1-dr))
	y
}

DE_ul <- function(C,re,b,d,dr,R_resi){
	y <- numeric()
	y[1] <- (factorial(2*re)/factorial(re))^dr*b*(re*(C+b-d))^(re*(1-dr))-(C-b-d)*R_resi^re
	y[2] <- d*R_resi^re-factorial(re)^dr*(C-b-d)*(re*(C+b-d))^(re*(1-dr))
	y
}

get.init_ul_cp <- function(C_seq,R,dr_seq,re_seq,death){
	d_seq <- C_seq*death
	b_seq <- 0.9*R*C_seq/crossprod(C_seq,re_seq)+d_seq-C_seq
	while(sum(b_seq<=0)|sum(C_seq-b_seq-d_seq)<=0&death<1){
		death <- death + 0.01
		d_seq <- C_seq*death
		b_seq <- 0.9*R*C_seq/crossprod(C_seq,re_seq)+d_seq-C_seq
	}
	c(b_seq,d_seq)
}

#solve_ul_cp(c(10,10),30,c(0.5,0.5),c(1,1))
solve_ul_cp <- function(C_seq,R,dr_seq,re_seq){##calculate growth and death rate given current abundances;set the initial death rate to be 0.1
	#re_seq <- re_seq/crossprod(C_seq,re_seq)*R  ##calculate actual resource requirement in the defined time frame
	death <- 0.01
	init.x <- get.init_ul_cp(C_seq,R,dr_seq,re_seq,death)
	solution <- nleqslv(init.x,de_series_ul_cp,C_seq=C_seq,R=R,dr_seq=dr_seq,re_seq=re_seq,control=list(maxit=1000))
	while(solution[[3]]!=1&death<1){
		#print(dr)
		death <- death+0.01
		init.x <- get.init_ul_cp(C_seq,R,dr_seq,re_seq,death)
		solution <- nleqslv(init.x,de_series_ul_cp,C_seq=C_seq,R=R,dr_seq=dr_seq,re_seq=re_seq,control=list(maxit=1000))
	}
	solved <- solution[[3]]==1|solution[[3]]==2
	if(!solved){
		print(solution[[4]])
		solution[[1]] <- rep(NA,4)## to get rid of the points that cannot be correctly solved
	}
	solution <- solution[[1]][1:length(C_seq)]-solution[[1]][1:length(C_seq)+length(C_seq)]
	solution <- c(solution,solution/C_seq)
	if(solved&out_of_bound(solution[-(1:length(C_seq))])){
			print("solution out of bound")
			solution <- rep(NA,2*length(C_seq))
		}
	solution
}

de_series_li_cp_beta<- function(x,C_seq,R,dr_seq,re_seq){## x[1:7] is b1/N1, d1/N1, b2/N2,..., E
	b_seq <- x[1:((length(x)-1)/2)]
	d_seq <- x[((length(x)-1)/2+1):(length(x)-1)]
	y <- numeric()
	for(i in 1:length(C_seq)){
		y[(2*i-1):(2*i)] <- DE(C_seq[i],re_seq[i],b_seq[i],d_seq[i],dr_seq[i],x[2*length(C_seq)+1])
	}
	Cprime_seq <- C_seq+b_seq-d_seq
	y[2*length(C_seq)+1] <- R-crossprod(Cprime_seq,re_seq)
	y
}


DE <- function(C,re,b,d,dr,E){##
	la <- -log(E)
	y <- numeric(2)
	y[1] <- (4*re/exp(1))^(re*dr)*b*(re*(C+b-d))^(re*(1-dr))-(C-b-d)*E^re
	y[2] <- d*(re/exp(1))^(dr*re)*E^re-(C-b-d)*(re*(C+b-d))^(re*(1-dr))
	y
}

DE <- function(C,re,b,d,dr,E){##corrected; E = exp(-la)
	y <- numeric(2)
	y[1] <- (factorial(2*re)/factorial(re))^dr*b*(re*(C+b-d))^(re*(1-dr))-(C-b-d)*E^re
	y[2] <- d*E^re-(C-b-d)*(re*(C+b-d))^(re*(1-dr))*factorial(re)^dr
	y
}

get.init_li_cp_beta <- function(C_seq,R,dr_seq,re_seq,death){
	d_seq <- C_seq*death
	b_seq <- R*C_seq/crossprod(C_seq,re_seq)+d_seq-C_seq
	while(sum(b_seq<0)|sum(C_seq-b_seq-d_seq)<=0&death<1){
		death <- death + 0.01
		d_seq <- C_seq*death
		b_seq <- R*C_seq/crossprod(C_seq,re_seq)+d_seq-C_seq
	}  
	E <- mean(b_seq*(C_seq+b_seq-d_seq)^(1-dr_seq)*2^dr_seq/(C_seq-b_seq-d_seq))
	if(E<0){
		E <- 1
	}
	c(b_seq,d_seq,E)
}

solve_li_cp_beta <- function(C_seq,R,dr_seq,re_seq){##calculate growth and death rate given current abundances;set the initial death rate to be 0.1
	if(R>2*(crossprod(C_seq,re_seq))){
		print("R too big")
		solution <- rep(NA,8)
	}else{
		death <- 0
		init.x <- get.init_li_cp_beta(C_seq,R,dr_seq,re_seq,death)
		solution <- nleqslv(init.x,de_series_li_cp_beta,C_seq=C_seq,R=R,dr_seq=dr_seq,re_seq=re_seq,control=list(xtol=1e-5,maxit=100,cndtol=1e-15))
		while(solution[[3]]!=1&death<1){
			death <- death+0.01
			init.x <- get.init_li_cp_beta(C_seq,R,dr_seq,re_seq,death)
		solution <- nleqslv(init.x,de_series_li_cp_beta,C_seq=C_seq,R=R,dr_seq=dr_seq,re_seq=re_seq)
		}
		solved <- solution[[3]]==1
		if(!solved){
			print(solution[[4]])
			solution[[1]] <- rep(NA,2*length(C_seq))## to get rid of the points that cannot be correctly solved
		}
		solution <- solution[[1]][1:length(C_seq)]-solution[[1]][1:length(C_seq)+length(C_seq)]
		solution <- c(solution,solution/C_seq)
		if(solved&out_of_bound(solution[-(1:length(C_seq))])){
			print("solution out of bound")
			solution <- rep(NA,2*length(C_seq))
		}
	}
		solution
}

out_of_bound <- function(sequence){ ## to tell whether rates solved are reasonable; not sure how to set
	sum(sequence< -1|sequence>1)>0
	#sum(sequence<0)>0
	#F
}

dynamics <- function(C_seq,R,dr_seq,re_seq,nsteps){
	time_track <- matrix(nrow=nsteps,ncol=length(C_seq))
	i <- 1
	while(i <=nsteps){
		time_track[i,] <- C_seq + solve_li_cp_beta(C_seq,R,dr_seq,re_seq)[1:length(C_seq)]
		C_seq <- time_track[i,]
		print(i)
		if(is.na(time_track[i,1])){
			i <- nsteps
		}else{
			i <- i+1
		}
	}
	time_track <- data.frame(time_track)
	names(time_track) <- paste("N",1:length(C_seq),sep="")
	time_track
}

#C_seq <- rep(10,3)
#R <- 6
#dr_seq <- rep(0.9,3)
#re_seq <- c(0.1,0.2,0.3)
#nmax <- 1000

#find_equi(C_seq,R,dr_seq,re_seq,nmax,0.001)

find_equi <- function(C_seq,R,Dr_seq,theta_seq,nmax=500,equi_judge=0.01){
	time_track <- matrix(nrow=nmax,ncol=length(C_seq))
	theta_seq <- theta_seq/crossprod(C_seq,theta_seq)*R#get the actual values of resource requirement
	i <- 1
	while(i <= nmax){
		increment <- try(solve_li_cp_beta(C_seq,R,Dr_seq,theta_seq)[1:length(C_seq)])
		if(class(increment)=="try-error"){
			break
		}
		time_track[i,] <- C_seq + increment
		C_seq <- time_track[i,]
		if(is.na(time_track[i,1])|max(abs(increment/C_seq))<equi_judge){##assume equilibrium is reached when maximum abundance change is smaller than 1%
			break
		}
		print(i)
		i <- i+1
	}
	print("one community finished")
	time_track[sum(!is.na(time_track[,1])),]
}


plot_dynamics <- function(time_track,ltyp){
	#limit<- c(min(time_track),max(time_track))
	time <- 1:dim(time_track)[1]-1
	#plot(time_track[,1]~time,xlab="time",ylab="abundance",ylim=limit)
	for(i in 1:dim(time_track)[2]){
		lines(time,time_track[,i],lty=ltyp)
	}
	#legend("bottomright",paste("Sp",1:dim(time_track)[2], " theta=",re_seq, " Dr= ",dr_seq,sep=""),pch=1,col=1:dim(time_track)[2])
}
