##predator can be partial or complete

predator_prey <- function(inits,theta,rs,SSp,nstep,partial=T){
	predator <- inits[1]
	prey <- inits[2]
	for(t in 1:nstep){
		##predator has to be complete redistribution
		dpredator <- rs[1]*predator[t]*(1-theta*(rs[1]+1)/rs[1]*predator[t]/prey[t])/(1+(rs[1]+1)*theta*predator[t]/prey[t])
		predator[t+1] <- predator[t]+ dpredator
		##but prey can be both
		if(partial){
			prey[t+1] <- prey[t]+ prey[t]*rs[2]*(1-prey[t]/SSp)/(1+(rs[2]-1)*prey[t]/SSp) - predator[t] * theta 
		}else{
			prey[t+1] <- prey[t]+ prey[t]*rs[2]*(1-prey[t]/SSp)/(1+rs[2]*prey[t]/SSp) - predator[t] * theta 
		}	
	}
	time <- 1: (nstep+1)
	plot(predator~time,col=2,ylim=range(c(predator,prey)))
	lines(time,predator,col=2)
	lines(time,prey,col=3)
	points(time,prey,col=3)
	legend("bottomright",c("predator","prey"),col=2:3,pch=1)
	list(predator= predator,prey = prey)
}

inits <- c(1,100)
theta <- 2
rs <- c(100,500)
SSp <- 1000
nstep <- 100
predator_prey(inits,theta,rs,SSp,nstep,partial=T)