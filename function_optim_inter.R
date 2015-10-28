##Get abundance, theta and CV from data

abundance <- function(sp_list,data){
	N_seq <- c()
	for(sp in sp_list){
		N_seq <- c(N_seq,sum(data$sp==sp&!is.na(data$dbh),na.rm=T))
	}
	N_seq
}

BA_CV <- function(sp_list,data,nsample){
	BA <- CV <- c() ##basal area and coefficient of variation
	for(sp in sp_list){
		dbh <- data$dbh[data$sp==sp&!is.na(data$dbh)]
		if(!is.na(nsample)){
			dbh <- sample(dbh,nsample)
		}
		ba <- dbh^2
		BA <- c(BA,mean(ba))
		CV <- c(CV,sd(dbh)/mean(dbh))
	}
	print(paste("correlation between CV and BA",cor.test(BA,CV)$p.value<0.05))
	list(BA=BA,CV=CV)
}

##function fitting

library(nleqslv)
Ns <- c(10,5,10,9)
BAs <- c(1,1.5,1,1.1)
CVs <- c(2,2,2,2.2)
inter <- T

Ns <- round(rnorm(30,100,10))
BAs <- 200-Ns+ rnorm(length(Ns),100,10)
CVs <- Ns/mean(Ns) +runif(length(Ns),1,2)

fit_MERA <- function(Ns,BAs,CVs,inter){
	R0 = sum(Ns*BAs)
	inits <- c(0.5,1)
	fit <- optim(inits,SSE,Ns=Ns,BAs=BAs,CVs=CVs,R0=R0,inter=inter)
	if(fit[[4]]>0){
		print(fit[[5]])
	}
	paras <- fit[[1]]
	MNs <- MERA_N(paras,BAs,CVs,R0,inter)
	plot(MNs~Ns,xlab="abundance (observed)",ylab="abundance (predicted)")
	abline(a=0,b=1,col=2)
	Y <- lm(MNs~Ns+1)
	R2 <- summary(Y)$r.squared
	legend("bottomright",c(paste(c("k=","mu="),round(paras,2)),paste("R2=",R2)))
	list(MNs=MNs,paras=paras)
}

SSE <- function(paras,Ns,BAs,CVs,R0,inter){
	MNs <- MERA_N(paras,BAs,CVs,R0,inter)
	sum((MNs-Ns)^2)
}

MERA_N <- function(paras,BAs,CVs,R0,inter=F){
	MPs <- transform(paras,BAs,CVs)
	Ts <- MPs$Ts
	Ds <- MPs$Ds
	Cs <- get_C(Ts,Ds,R0,inter)
	Cs^Ds*Ts/BAs
}

transform <- function(paras,BAs,CVs){
	Drs <- paras[1]*CVs + paras[2]
	Ds <- 1/(Drs-1)
	Ts <- BAs^(Drs*Ds)
	list(Ts=Ts,Ds=Ds)
}

transform <- function(paras,BAs,CVs){
	Drs <- logi(CVs,paras[1],paras[2])
	Ds <- 1/(Drs-1)
	Ts <- BAs^(Drs*Ds)
	list(Ts=Ts,Ds=Ds)
}

logi <- function(x,k,mu){
	1/(1+exp(-k*(x-mu)))
}

get_C <- function(Ts,Ds,R0,inter){
	if(inter){
		Cs.init <- exp(log(R0/sum(Ts))/Ds)
	}else{
		Cs.init <- exp(log(R0/sum(Ts))/mean(Ds))
	}
	#print(Cs.init)
	Cs <- nleqslv(Cs.init,constraints,Ts=Ts,Ds=Ds,R0=R0,inter=inter)
	if(Cs[[3]]!=1){
		print(Cs[[4]])
		print(Cs[[2]])
	}
	Cs[[1]]
}

constraints <- function(Cs,Ts,Ds,R0,inter){
	Y <- c()
	if(!inter){
		Y = sum(Cs^Ds*Ts)-R0
	}else{
		for(i in 1:(length(Cs)-1)){
			Y <- c(Y,(sum(Cs[-i]^Ds[i]*Ts[i]) - sum(Cs[i]^Ds[-i]*Ts[-i])))
		} 
		Y = c(Y,(sum(Cs^Ds*Ts)-R0))
	}
	Y	
}